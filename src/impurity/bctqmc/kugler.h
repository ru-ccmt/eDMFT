#ifndef EDMFT_KUGLER_H
#define EDMFT_KUGLER_H
#include <array>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>
namespace edmft::kugler {
// Unnormalized node correlators. Each pair separates the two existing Nambu
// endpoint sectors: creation replacement, annihilation replacement, both.
// Empty storage is intentional when the optional estimator is disabled.
struct Correlators {
  std::array<Eigen::MatrixXd,6> value;
  bool empty() const { return value[0].size()==0; }
  void ensure(int rows,int cols) {
    if (empty()) for (auto& x:value) x=Eigen::MatrixXd::Zero(rows,cols);
    else if (value[0].rows()!=rows || value[0].cols()!=cols)
      throw std::logic_error("Kugler correlator dimensions changed");
  }
  void add(const Correlators& other,double scale=1.) {
    if (other.empty()) return;
    ensure(other.value[0].rows(),other.value[0].cols());
    for (int i=0;i<6;++i) value[i]+=scale*other.value[i];
  }
  void clear() { for (auto& x:value) x.resize(0,0); }
};
// Reconstruct physical orbital/Nambu matrices from the two legacy endpoint
// sectors. Particle-hole conjugation swaps LEFT and RIGHT mixed correlators.
// Hartree is an equal-time moment: its hole sign differs from the frequency
// correlator, since conjugation also reverses the leading 1/(iw) term.
template<class CixType,class Blocks,class Value>
Eigen::MatrixXcd reconstruct(const CixType& cix,const Blocks& bb,int block,
    int kind,bool hartree,const Value& value) {
  const int n=bb.blocks[block].size();
  const bool sc=bb.is_block_sc[block];
  Eigen::MatrixXcd out=Eigen::MatrixXcd::Zero(n,n);
  for (int b=0;b<cix.Nbath;++b) if (bb.orb_ind[b][0]==block) {
    int row=bb.orb_ind[b][1],col=bb.orb_ind[b][2];
    const int s0=cix.bath_sigind0(b),s=std::abs(s0);
    const bool rh=sc && row>=n/2,ch=sc && col>=n/2;
    const int sector=ch ? 1 : 0;
    const int source_kind=rh && ch && kind<2 && !hartree ? 1-kind : kind;
    std::complex<double> v=value(2*source_kind+sector,s);
    v*=s0<0 ? -1. : 1.;
    if (rh && ch) { v=-std::conj(v); std::swap(row,col); }
    if ((kind==0 || kind==2) && ch) v=-v;
    if ((kind==1 || kind==2) && rh) v=-v;
    if (hartree && rh && ch) v=-v;
    out(row,col)=v;
  }
  return out;
}
inline Eigen::MatrixXcd self_energy(const Eigen::MatrixXcd& hartree,
    const Eigen::MatrixXcd& interaction,const Eigen::MatrixXcd& left,
    const Eigen::MatrixXcd& green,const Eigen::MatrixXcd& right) {
  // Solve rather than explicitly invert G. Order matters for orbital/Nambu matrices.
  Eigen::MatrixXcd result=hartree+interaction-left*green.partialPivLu().solve(right);
  if (!result.allFinite()) throw std::runtime_error("nonfinite Kugler self-energy");
  return result;
}
}
#endif
