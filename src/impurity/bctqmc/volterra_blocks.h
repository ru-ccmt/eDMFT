#ifndef BCTQMC_VOLTERRA_BLOCKS_H
#define BCTQMC_VOLTERRA_BLOCKS_H

#include <numeric>
#include <vector>

namespace volterra_blocks {
using Parts=std::vector<std::vector<int>>;

// Use structural zeros, not values of Sigma at a particular iteration. A
// nonzero Hamiltonian entry is an edge even when Sindex does not store it.
inline Parts partition(const Array<double,2>& h, const Array<int,2>& index,
                        bool enabled=true)
{
  const int n=h.extent(0);
  std::vector<int> parent(n); std::iota(parent.begin(),parent.end(),0);
  auto root=[&](int i) {
    while(parent[i]!=i) { parent[i]=parent[parent[i]]; i=parent[i]; }
    return i;
  };
  if (!enabled) { Parts all(1); all[0]=parent; return all; }
  for (int i=0; i<n; ++i) for (int j=0; j<i; ++j)
    if (h(i,j)!=0.0 || h(j,i)!=0.0 || index(i,j)!=0 || index(j,i)!=0) {
      const int a=root(i),b=root(j); parent[std::max(a,b)]=std::min(a,b);
    }
  Parts parts;
  std::vector<int> position(n,-1);
  for (int i=0; i<n; ++i) {
    const int r=root(i);
    if (position[r]<0) { position[r]=parts.size(); parts.emplace_back(); }
    parts[position[r]].push_back(i);
  }
  return parts;
}

struct Result {
  std::vector<double> probability;
  double diagonalization=0.0, sigma_rotation=0.0, volterra=0.0,
         g_rotation=0.0, scatter=0.0;
};

// Solve independent pieces in the original product basis, with the unchanged
// scalar/matrix Volterra kernels supplied by the caller. Emit only one piece
// at a time; no full-container time-dependent matrix is allocated here.
template<class ScalarSolve,class MatrixSolve,class Emit>
Result solve(const Array<double,2>& ham, const Array<int,2>& index,
             const Array<double,2>& sigma, const Array<double,1>& tau,
             double lambda, const Parts& parts,
             ScalarSolve scalar_solve, MatrixSolve matrix_solve, Emit emit)
{
  Result result;
  const int n=ham.extent(0), nt=tau.size();
  Array<double,2> endpoint(n,n); endpoint=0.0;
  for (const auto& part : parts) {
    const int d=part.size();
    Array<double,3> g;
    if (d==1) {
      const int i=part[0],label=index(i,i);
      if (!label) throw std::runtime_error("Volterra singleton has zero diagonal Sindex");
      Array<double,1> values;
      const double t0=get_cpu_time();
      Array<double,1> self(sigma(std::abs(label)-1,Range::all()).copy());
      if (label<0) self*=-1;
      scalar_solve(values,self,ham(i,i)-lambda);
      result.volterra+=get_cpu_time()-t0;
      g.resize(nt,1,1); g(Range::all(),0,0)=values;
    } else {
      const double t0=get_cpu_time();
      Array<double,2> h(d,d),s1(d,d),s2(d,d);
      for (int i=0; i<d; ++i) for (int j=0; j<d; ++j) h(i,j)=ham(part[i],part[j]);
      Array<double,1> e(d),work(3*d);
      if (xsyev(d,h.data(),d,e.data(),work.data(),work.size(),"V"))
        throw std::runtime_error("Volterra subblock diagonalization failed");
      const double t1=get_cpu_time(); result.diagonalization+=t1-t0;
      Array<double,3> sg(d,d,nt);
      for (int it=0; it<nt; ++it) {
        s1=0.0;
        for (int i=0; i<d; ++i) for (int j=0; j<d; ++j) {
          const int label=index(part[i],part[j]);
          if (label) s1(i,j)=sign(label)*sigma(std::abs(label)-1,it);
        }
        tiny::dot(s2,h,s1);
        tiny::dot(s1,s2,h,false,true);
        sg(Range::all(),Range::all(),it)=s1;
      }
      const double t2=get_cpu_time(); result.sigma_rotation+=t2-t1;
      e-=lambda;
      matrix_solve(g,sg,e);
      const double t3=get_cpu_time(); result.volterra+=t3-t2;
      for (int it=0; it<nt; ++it) {
        tiny::dot(s1,h,g(it,Range::all(),Range::all()),true,false);
        tiny::dot(s2,s1,h);
        g(it,Range::all(),Range::all())=s2;
      }
      result.g_rotation+=get_cpu_time()-t3;
    }
    for (int i=0; i<d; ++i) for (int j=0; j<d; ++j)
      endpoint(part[i],part[j])=g(nt-1,i,j);
    const double t0=get_cpu_time(); emit(part,g); result.scatter+=get_cpu_time()-t0;
  }

  // Preserve the original full-Hamiltonian eigenstate ordering and convention
  // of Probability.dat, including degenerate states in different pieces. The
  // full eigenbasis is needed only at this endpoint, never for time stepping.
  const double t0=get_cpu_time();
  Array<double,2> full(ham.copy());
  Array<double,1> energies(n),work(3*n);
  if (xsyev(n,full.data(),n,energies.data(),work.data(),work.size(),"V"))
    throw std::runtime_error("Volterra probability diagonalization failed");
  result.diagonalization+=get_cpu_time()-t0;
  result.probability.assign(n,0.0);
  for (int k=0; k<n; ++k)
    for (const auto& part : parts)
      for (int i : part) {
        double value=0.0;
        for (int j : part) value+=endpoint(i,j)*full(k,j);
        result.probability[k]+=full(k,i)*value;
      }
  return result;
}
} // namespace volterra_blocks
#endif
