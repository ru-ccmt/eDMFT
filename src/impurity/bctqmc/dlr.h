// Lightweight fermionic Discrete Lehmann Representation (DLR).
//
// The construction follows the interpolative-decomposition algorithm of
// J. Kaye, K. Chen, and O. Parcollet, Phys. Rev. B 105, 235115 (2022),
// using composite Gaussian/Chebyshev grids and pivoted reorthogonalized
// Gram-Schmidt.  The implementation is self-contained apart from Eigen.
#ifndef EDMFT_DLR_H
#define EDMFT_DLR_H

#include "Eigen/Dense"

#include <algorithm>
#include <cmath>
#include <complex>
#include <chrono>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace edmft::dlr {

using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::VectorXd;

// Standard fermionic continuation kernel, with dimensionless beta=1:
// K(t,w) = -exp(-t*w)/(1+exp(-w)), 0 <= t <= 1.
inline double kernel_abs(double t, double w)
{
  if (w >= 0.0) return -std::exp(-t*w)/(1.0+std::exp(-w));
  return -std::exp((1.0-t)*w)/(1.0+std::exp(w));
}

// Relative time keeps full relative precision close to both endpoints:
// t in [0,1/2] represents itself and t in [-1/2,0) represents 1+t.
inline double kernel_relative(double t, double w)
{
  return t >= 0.0 ? kernel_abs(t,w) : kernel_abs(-t,-w);
}

inline double absolute_to_relative(double t)
{
  return (t > 0.5 && t < 1.0) ? t-1.0 : t;
}

inline double relative_to_absolute(double t)
{
  return t < 0.0 ? t+1.0 : t;
}

struct Options {
  // Physical inverse temperature.  Internally the basis is constructed with
  // dimensionless t=tau/beta and w=beta*omega.
  double beta = 1.0;
  // The pseudo spectral window.  lambda=beta*omega_max controls the numerical
  // rank much more directly than beta and omega_max separately.
  double omega_max = 1.0;
  // Requested kernel compression tolerance, not a guaranteed error for an
  // arbitrary input function whose spectrum extends outside omega_max.
  double eps = 1e-10;
  // Nodes per dyadic panel on the temporary fine grid.  This controls how
  // accurately the continuous kernel space is sampled before compression;
  // it is not the final DLR rank.
  int panel_order = 24;
};

namespace detail {

inline std::pair<VectorXd,VectorXd> gauss_legendre(int n)
{
  VectorXd x(n), w(n);
  const int nhalf = (n+1)/2;
  for (int i=0; i<nhalf; ++i) {
    double z = std::cos(std::numbers::pi*(i+0.75)/(n+0.5));
    for (int it=0; it<100; ++it) {
      double p0=1.0, p1=z;
      for (int k=2; k<=n; ++k) {
        const double p2=((2.0*k-1.0)*z*p1-(k-1.0)*p0)/k;
        p0=p1; p1=p2;
      }
      const double pn = n==0 ? 1.0 : (n==1 ? z : p1);
      const double pnm1 = n<=1 ? 1.0 : p0;
      const double dp = n*(z*pn-pnm1)/(z*z-1.0);
      const double dz = -pn/dp;
      z += dz;
      if (std::abs(dz)<2e-15) break;
    }
    double p0=1.0, p1=z;
    for (int k=2; k<=n; ++k) {
      const double p2=((2.0*k-1.0)*z*p1-(k-1.0)*p0)/k;
      p0=p1; p1=p2;
    }
    const double pn = n==1 ? z : p1;
    const double pnm1 = n<=1 ? 1.0 : p0;
    const double dp = n*(z*pn-pnm1)/(z*z-1.0);
    const double wi = 2.0/((1.0-z*z)*dp*dp);
    x(i)=-z; x(n-1-i)=z;
    w(i)=wi; w(n-1-i)=wi;
  }
  return {x,w};
}

struct FineGrid {
  int p;
  int frequency_panels;
  int time_panels;
  VectorXd omega;
  VectorXd time_relative;
  VectorXd sqrt_time_weights;
};

inline FineGrid make_fine_grid(double lambda, int p)
{
  if (!(lambda>0.0) || p<4) throw std::invalid_argument("DLR requires lambda>0 and panel_order>=4");
  // Dyadic panels resolve the sharp endpoint structure that develops at
  // large lambda.  The final number of retained frequencies is selected by
  // numerical rank and is normally far smaller than this candidate grid.
  const int npom=std::max(static_cast<int>(std::ceil(std::log2(lambda))),1);
  const int npt=std::max(static_cast<int>(std::ceil(std::log2(lambda)))-2,1);

  VectorXd cheb(p);
  for (int i=0; i<p; ++i)
    cheb(i)=0.5*(1.0-std::cos((2.0*i+1.0)*std::numbers::pi/(2.0*p)));

  VectorXd omega(2*p*npom);
  double a=0.0;
  for (int panel=0; panel<npom; ++panel) {
    const double b=lambda/std::pow(2.0,npom-panel-1);
    for (int i=0; i<p; ++i) omega((npom+panel)*p+i)=a+(b-a)*cheb(i);
    a=b;
  }
  for (int i=0; i<p*npom; ++i) omega(i)=-omega(2*p*npom-1-i);

  auto [xgl,wgl]=gauss_legendre(p);
  xgl=(xgl.array()+1.0)/2.0;
  VectorXd time(2*p*npt), sqrtw(2*p*npt);
  a=0.0;
  for (int panel=0; panel<npt; ++panel) {
    const double b=1.0/std::pow(2.0,npt-panel);
    for (int i=0; i<p; ++i) {
      time(panel*p+i)=a+(b-a)*xgl(i);
      sqrtw(panel*p+i)=std::sqrt(0.5*(b-a)*wgl(i));
    }
    a=b;
  }
  for (int i=0; i<p*npt; ++i) {
    time(p*npt+i)=-time(p*npt-1-i);
    sqrtw(p*npt+i)=sqrtw(p*npt-1-i);
  }
  return {p,npom,npt,std::move(omega),std::move(time),std::move(sqrtw)};
}

// Select rows by pivoted, reorthogonalized Gram-Schmidt.  If target_rank is
// positive, select exactly that many rows; otherwise stop at residual eps.
inline std::vector<int> select_rows(MatrixXd a, double eps, int target_rank=-1)
{
  const int m=a.rows(), n=a.cols();
  const int maxrank=target_rank>0 ? std::min({target_rank,m,n}) : std::min(m,n);
  VectorXd norms(m);
  for (int i=0; i<m; ++i) norms(i)=a.row(i).squaredNorm();
  std::vector<int> piv(m);
  for (int i=0; i<m; ++i) piv[i]=i;
  std::vector<int> selected;
  selected.reserve(maxrank);
  const double eps2=eps*eps;

  for (int j=0; j<maxrank; ++j) {
    int jpiv=j;
    for (int k=j+1; k<m; ++k) if (norms(k)>norms(jpiv)) jpiv=k;
    if (jpiv!=j) {
      a.row(j).swap(a.row(jpiv));
      std::swap(norms(j),norms(jpiv));
      std::swap(piv[j],piv[jpiv]);
    }
    for (int k=0; k<j; ++k) a.row(j)-=a.row(j).dot(a.row(k))*a.row(k);
    const double nrm=a.row(j).squaredNorm();
    if (target_rank<0 && nrm<=eps2) break;
    if (!(nrm>64.0*std::numeric_limits<double>::min()))
      throw std::runtime_error("DLR node selection lost numerical rank");
    a.row(j)/=std::sqrt(nrm);
    selected.push_back(piv[j]);
    for (int k=j+1; k<m; ++k) {
      if (target_rank<0 && norms(k)<=eps2) continue;
      a.row(k)-=a.row(k).dot(a.row(j))*a.row(j);
      norms(k)=a.row(k).squaredNorm();
    }
  }
  return selected;
}

} // namespace detail

class FermionicBasis {
public:
  explicit FermionicBasis(const Options& opts)
    : beta_(opts.beta), omega_max_(opts.omega_max), eps_(opts.eps),
      lambda_(opts.beta*opts.omega_max)
  {
    if (!(beta_>0.0) || !(omega_max_>0.0) || !(eps_>0.0))
      throw std::invalid_argument("DLR requires beta, omega_max, and eps to be positive");
    build(opts.panel_order);
  }

  // Number r of selected frequencies and selected time nodes.  All compact
  // transforms below are r-by-r; the original dense tau mesh is used only to
  // fit input data or to evaluate final output.
  int rank() const { return static_cast<int>(omega_dimensionless_.size()); }
  double beta() const { return beta_; }
  double omega_max() const { return omega_max_; }
  double lambda() const { return lambda_; }
  double tolerance() const { return eps_; }
  // Conservative gate for the optional coefficient-chain optimization.
  // A small residual is necessary, not a general conditioning guarantee.
  // Tight interpolation tolerances were less stable in the SC stress tests;
  // retain the node/solve route there even if the caller requests chaining.
  bool coefficient_convolution_available() const { return coefficient_convolution_ok_; }
  const VectorXd& dimensionless_frequencies() const { return omega_dimensionless_; }
  const VectorXd& frequencies() const { return omega_physical_; }
  const VectorXd& relative_time_nodes() const { return time_relative_; }
  const VectorXd& time_nodes() const { return time_physical_; }
  const MatrixXd& coefficients_to_values_matrix() const { return cf2it_; }
  const MatrixXd& reflection_values_matrix() const { return reflection_values_; }

  // Reflect values sampled at the representative nodes from tau_k to
  // beta-tau_k.  The dense r-by-r map is constructed once with the basis;
  // diagram evaluators therefore need neither a new LU solve nor r separate
  // kernel evaluations for every pair family and every AAA pole.
  MatrixXd reflect_values(const MatrixXd& values) const
  {
    check_rows(values.rows());
    return reflection_values_*values;
  }

  // Pair contractions naturally store one pair per row and one DLR time node
  // per column.  In that orientation reflection is a right multiplication.
  // This form is especially useful for future MPI row-block decomposition:
  // every rank can reflect its owned pair rows independently.
  MatrixXd reflect_transposed_values(const MatrixXd& transposed_values) const
  {
    check_rows(transposed_values.cols());
    return transposed_values*reflection_values_.transpose();
  }

  VectorXd values_to_coefficients(const VectorXd& values) const
  {
    check_rows(values.size());
    return lu_.solve(values);
  }
  VectorXcd values_to_coefficients(const VectorXcd& values) const
  {
    check_rows(values.size());
    return zlu_.solve(values);
  }
  MatrixXd values_to_coefficients(const MatrixXd& values) const
  {
    check_rows(values.rows());
    // cf2it_ is square by construction: representative time nodes and DLR
    // frequencies have the same rank.  The LU factorization is built once.
    return lu_.solve(values);
  }
  MatrixXcd values_to_coefficients(const MatrixXcd& values) const
  {
    check_rows(values.rows());
    return zlu_.solve(values);
  }

  VectorXd coefficients_to_values(const VectorXd& coeffs) const
  {
    check_rows(coeffs.size());
    return cf2it_*coeffs;
  }
  VectorXcd coefficients_to_values(const VectorXcd& coeffs) const
  {
    check_rows(coeffs.size());
    return cf2it_*coeffs;
  }
  MatrixXd coefficients_to_values(const MatrixXd& coeffs) const
  {
    check_rows(coeffs.rows());
    return cf2it_*coeffs;
  }
  MatrixXcd coefficients_to_values(const MatrixXcd& coeffs) const
  {
    check_rows(coeffs.rows());
    return cf2it_*coeffs;
  }

  // Apply the divided-difference/Hilbert matrix used by the analytic
  // convolution formula to one or many coefficient vectors.  Exposing this
  // small operation lets a pole-major diagram evaluator prepare H*g once
  // when the same right-hand propagator g is reused for many bath pairs.
  MatrixXd convolution_hilbert(const MatrixXd& coeffs) const
  {
    check_rows(coeffs.rows());
    return thilb_*coeffs;
  }

  template <typename Scalar>
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>
  fit_values(const VectorXd& tau, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& values) const
  {
    if (tau.size()!=values.rows()) throw std::invalid_argument("DLR fit: tau/value size mismatch");
    MatrixXd k(tau.size(),rank());
    for (int i=0; i<tau.size(); ++i) {
      if (tau(i)<0.0 || tau(i)>beta_) throw std::out_of_range("DLR fit: tau outside [0,beta]");
      for (int l=0; l<rank(); ++l) k(i,l)=kernel_abs(tau(i)/beta_,omega_dimensionless_(l));
    }
    // This is the only intentionally overdetermined transform: it projects a
    // function sampled on the solver's dense/nonuniform tau mesh into the r
    // DLR coefficients.  Subsequent transforms use the square node matrix.
    if constexpr (std::is_same_v<Scalar,std::complex<double>>)
      return k.cast<std::complex<double> >().colPivHouseholderQr().solve(values);
    else
      return k.colPivHouseholderQr().solve(values);
  }

  template <typename Scalar>
  Eigen::Matrix<Scalar,1,Eigen::Dynamic>
  evaluate(double tau, const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>& coeffs) const
  {
    check_rows(coeffs.rows());
    if (tau<0.0 || tau>beta_) throw std::out_of_range("DLR evaluate: tau outside [0,beta]");
    Eigen::RowVectorXd k(rank());
    for (int l=0; l<rank(); ++l) k(l)=kernel_abs(tau/beta_,omega_dimensionless_(l));
    return k*coeffs;
  }

  // Values at DLR nodes of
  //
  //        h(tau) = int_0^tau f(tau-t') g(t') dt'.
  //
  // Inputs are DLR coefficients, not time-node values.  For unequal DLR
  // frequencies the analytic integral contains a divided difference; those
  // terms are encoded by thilb_.  Equal frequencies give tau*K and are
  // encoded separately by ttcf2it_.  This is why the result remains in the
  // same r-dimensional representation without introducing tau*K as a second
  // permanent basis.  Cost is O(r^2) for one scalar convolution.
  VectorXcd time_ordered_convolution_values(const VectorXcd& f,
                                             const VectorXcd& g) const
  {
    check_rows(f.size()); check_rows(g.size());
    const VectorXcd diag=(f.array()*g.array()).matrix();
    const VectorXcd hf=thilb_.cast<std::complex<double> >()*f;
    const VectorXcd hg=thilb_.cast<std::complex<double> >()*g;
    const VectorXcd off=(f.array()*hg.array()+g.array()*hf.array()).matrix();
    return beta_*(ttcf2it_.cast<std::complex<double> >()*diag+
                  cf2it_.cast<std::complex<double> >()*off);
  }

  VectorXcd time_ordered_convolution_coefficients(const VectorXcd& f,
                                                   const VectorXcd& g) const
  {
    return values_to_coefficients(time_ordered_convolution_values(f,g));
  }

  // Batched real-valued counterpart of time_ordered_convolution_values.
  // Each column of g contains one set of DLR coefficients.  The left-hand
  // function f is shared by all columns, which lets Eigen evaluate the
  // Hilbert transforms and the final DLR transforms as dense matrix-matrix
  // products.  This is the hot operation in compact diagram evaluation.
  //
  // Algebraically, for every column j:
  //   diag_j = f .* g_j
  //   off_j  = f .* (H g_j) + (H f) .* g_j
  //   h_j    = beta [T diag_j + K off_j].
  // H=thilb_, T=ttcf2it_, and K=cf2it_.  Writing it in this order shares Hf
  // and turns r*Npair scalar operations into level-3 Eigen matrix products.
  MatrixXd time_ordered_convolution_values(const VectorXd& f,
                                            const MatrixXd& g) const
  {
    check_rows(f.size()); check_rows(g.rows());
    const VectorXd hf=thilb_*f;
    const MatrixXd hg=thilb_*g;
    MatrixXd diag=g;
    MatrixXd off=hg;
    MatrixXd hfg=g;
    diag.array().colwise() *= f.array();
    off.array().colwise() *= f.array();
    hfg.array().colwise() *= hf.array();
    off += hfg;
    return beta_*(ttcf2it_*diag+cf2it_*off);
  }

  MatrixXd time_ordered_convolution_coefficients(const VectorXd& f,
                                                  const MatrixXd& g) const
  {
    return values_to_coefficients(time_ordered_convolution_values(f,g));
  }

  // Matrix-valued counterpart. f[l] and g[l] are coefficient matrices and
  // ordinary matrix multiplication is used under the convolution integral.
  std::vector<MatrixXcd>
  time_ordered_convolution_values(const std::vector<MatrixXcd>& f,
                                  const std::vector<MatrixXcd>& g) const
  {
    if (static_cast<int>(f.size())!=rank() || static_cast<int>(g.size())!=rank())
      throw std::invalid_argument("DLR matrix convolution: coefficient count != rank");
    if (f.empty() || f[0].cols()!=g[0].rows())
      throw std::invalid_argument("DLR matrix convolution: incompatible matrices");
    const int nr=f[0].rows(), nc=g[0].cols();
    std::vector<MatrixXcd> hf(rank(),MatrixXcd::Zero(f[0].rows(),f[0].cols()));
    std::vector<MatrixXcd> hg(rank(),MatrixXcd::Zero(g[0].rows(),g[0].cols()));
    for (int l=0; l<rank(); ++l)
      for (int k=0; k<rank(); ++k) {
        hf[l].noalias() += thilb_(l,k)*f[k];
        hg[l].noalias() += thilb_(l,k)*g[k];
      }
    std::vector<MatrixXcd> d(rank()), off(rank());
    for (int l=0; l<rank(); ++l) {
      d[l]=f[l]*g[l];
      off[l]=hf[l]*g[l]+f[l]*hg[l];
    }
    std::vector<MatrixXcd> out(rank(),MatrixXcd::Zero(nr,nc));
    for (int i=0; i<rank(); ++i)
      for (int l=0; l<rank(); ++l)
        out[i].noalias() += beta_*(ttcf2it_(i,l)*d[l]+cf2it_(i,l)*off[l]);
    return out;
  }

  // Real-valued matrix counterpart used by the TCA block-matrix evaluator.
  // The order of matrix multiplication is the order of the arguments: this
  // represents int_0^t f(t') g(t-t') dt'.  For scalar functions this is the
  // same convolution as above after t' -> t-t'; for matrices preserving the
  // product order is essential because f and g need not commute.
  std::vector<MatrixXd>
  time_ordered_convolution_values(const std::vector<MatrixXd>& f,
                                  const std::vector<MatrixXd>& g) const
  {
    if (static_cast<int>(f.size())!=rank() ||
        static_cast<int>(g.size())!=rank())
      throw std::invalid_argument(
        "DLR real-matrix convolution: coefficient count != rank");
    if (f.empty() || f[0].cols()!=g[0].rows())
      throw std::invalid_argument(
        "DLR real-matrix convolution: incompatible matrices");
    const int nr=f[0].rows(),nk=f[0].cols(),nc=g[0].cols();
    for (int l=0;l<rank();++l)
      if (f[l].rows()!=nr || f[l].cols()!=f[0].cols() ||
          g[l].rows()!=g[0].rows() || g[l].cols()!=nc)
        throw std::invalid_argument(
          "DLR real-matrix convolution: matrix dimensions vary with coefficient");

    // Flatten matrix entries into multiple right-hand sides.  The previous
    // implementation evaluated Hf, Hg, and the final DLR transforms as r^2
    // separate MatrixXd additions.  For the rank~50 matrix-valued TCA
    // convolutions that dispatch overhead dominated the useful work.  These
    // are ordinary linear transforms in the DLR index, so four level-3 dense
    // products perform exactly the same algebra much more efficiently.
    MatrixXd f_flat(rank(),nr*nk),g_flat(rank(),nk*nc);
    for (int l=0;l<rank();++l) {
      for (int i=0;i<nr;++i)
        for (int j=0;j<nk;++j)
          f_flat(l,i*nk+j)=f[l](i,j);
      for (int i=0;i<nk;++i)
        for (int j=0;j<nc;++j)
          g_flat(l,i*nc+j)=g[l](i,j);
    }
    MatrixXd hf_flat(rank(),nr*nk),hg_flat(rank(),nk*nc);
    hf_flat.noalias()=thilb_*f_flat;
    hg_flat.noalias()=thilb_*g_flat;

    MatrixXd diagonal_flat(rank(),nr*nc),off_diagonal_flat(rank(),nr*nc);
    for (int l=0;l<rank();++l) {
      MatrixXd hf_l(nr,nk),hg_l(nk,nc);
      for (int i=0;i<nr;++i)
        for (int j=0;j<nk;++j)
          hf_l(i,j)=hf_flat(l,i*nk+j);
      for (int i=0;i<nk;++i)
        for (int j=0;j<nc;++j)
          hg_l(i,j)=hg_flat(l,i*nc+j);
      const MatrixXd diagonal=f[l]*g[l];
      const MatrixXd off_diagonal=hf_l*g[l]+f[l]*hg_l;
      for (int i=0;i<nr;++i)
        for (int j=0;j<nc;++j) {
          diagonal_flat(l,i*nc+j)=diagonal(i,j);
          off_diagonal_flat(l,i*nc+j)=off_diagonal(i,j);
        }
    }
    MatrixXd values_flat(rank(),nr*nc);
    values_flat.noalias()=beta_*(ttcf2it_*diagonal_flat+
                                 cf2it_*off_diagonal_flat);
    std::vector<MatrixXd> values(
      rank(),MatrixXd::Zero(nr,nc));
    for (int l=0;l<rank();++l)
      for (int i=0;i<nr;++i)
        for (int j=0;j<nc;++j)
          values[l](i,j)=values_flat(l,i*nc+j);
    return values;
  }

  // The same real matrix convolution with H*g supplied by the caller.
  // Pole-major strong-coupling diagrams reuse a dressed pseudo propagator g
  // for many bath pairs.  Recomputing thilb_*g for every pair discarded much
  // of that reuse.  The left function f remains bath dependent, so H*f is
  // still formed here.  g and hg are coefficient-space matrix functions.
  std::vector<MatrixXd>
  time_ordered_convolution_values_right_prepared(
      const std::vector<MatrixXd>& f,
      const std::vector<MatrixXd>& g,
      const std::vector<MatrixXd>& hg,
      bool reflect_output=false) const
  {
    if (static_cast<int>(f.size())!=rank() ||
        static_cast<int>(g.size())!=rank() ||
        static_cast<int>(hg.size())!=rank())
      throw std::invalid_argument(
        "DLR prepared real-matrix convolution: coefficient count != rank");
    if (f.empty() || f[0].cols()!=g[0].rows() ||
        hg[0].rows()!=g[0].rows() || hg[0].cols()!=g[0].cols())
      throw std::invalid_argument(
        "DLR prepared real-matrix convolution: incompatible matrices");
    const int nr=f[0].rows(),nk=f[0].cols(),nc=g[0].cols();
    for (int l=0;l<rank();++l)
      if (f[l].rows()!=nr || f[l].cols()!=nk ||
          g[l].rows()!=nk || g[l].cols()!=nc ||
          hg[l].rows()!=nk || hg[l].cols()!=nc)
        throw std::invalid_argument(
          "DLR prepared real-matrix convolution: matrix dimensions vary");

    // Only f must be flattened for a DLR-index transform.  Both g and H*g
    // have already been prepared and remain in their natural block-matrix
    // representation for the noncommuting atomic products below.
    MatrixXd f_flat(rank(),nr*nk);
    for (int l=0;l<rank();++l)
      for (int i=0;i<nr;++i)
        for (int j=0;j<nk;++j)
          f_flat(l,i*nk+j)=f[l](i,j);
    MatrixXd hf_flat(rank(),nr*nk);
    hf_flat.noalias()=thilb_*f_flat;

    MatrixXd diagonal_flat(rank(),nr*nc),off_diagonal_flat(rank(),nr*nc);
    for (int l=0;l<rank();++l) {
      MatrixXd hf_l(nr,nk);
      for (int i=0;i<nr;++i)
        for (int j=0;j<nk;++j)
          hf_l(i,j)=hf_flat(l,i*nk+j);
      const MatrixXd diagonal=f[l]*g[l];
      const MatrixXd off_diagonal=hf_l*g[l]+f[l]*hg[l];
      for (int i=0;i<nr;++i)
        for (int j=0;j<nc;++j) {
          diagonal_flat(l,i*nc+j)=diagonal(i,j);
          off_diagonal_flat(l,i*nc+j)=off_diagonal(i,j);
        }
    }
    MatrixXd values_flat(rank(),nr*nc);
    if (reflect_output)
      values_flat.noalias()=beta_*(
        reflected_ttcf2it_*diagonal_flat+
        reflected_cf2it_*off_diagonal_flat);
    else
      values_flat.noalias()=beta_*(ttcf2it_*diagonal_flat+
                                   cf2it_*off_diagonal_flat);
    std::vector<MatrixXd> values(rank(),MatrixXd::Zero(nr,nc));
    for (int l=0;l<rank();++l)
      for (int i=0;i<nr;++i)
        for (int j=0;j<nc;++j)
          values[l](i,j)=values_flat(l,i*nc+j);
    return values;
  }

  // Batched variant of the prepared right convolution.  Every item may have
  // different atomic-matrix dimensions, but all of them share the same DLR
  // grid.  Concatenating their scalar matrix entries into wide rank-by-M
  // matrices lets the expensive DLR transforms run as one multiple-RHS
  // operation instead of one small solve/GEMM per symmetry-block pair.
  //
  // f_values is deliberately supplied at the interpolation nodes.  The
  // vertex immediately before a convolution is pointwise in tau, so this
  // routine folds the required values-to-coefficients solve into the batch.
  // g and hg are the already prepared coefficient-space pseudo propagator
  // and its convolution Hilbert transform.  The atomic products themselves
  // remain separate because their dimensions and multiplication order vary.
  std::vector<std::vector<MatrixXd>>
  time_ordered_convolution_values_right_prepared_batch_reference(
      const std::vector<const std::vector<MatrixXd>*>& f_values,
      const std::vector<const std::vector<MatrixXd>*>& g,
      const std::vector<const std::vector<MatrixXd>*>& hg,
      bool reflect_output=false) const
  {
    const int count=static_cast<int>(f_values.size());
    if (static_cast<int>(g.size())!=count ||
        static_cast<int>(hg.size())!=count)
      throw std::invalid_argument(
        "DLR batched prepared convolution: batch sizes differ");
    if (count==0) return {};

    std::vector<int> nr(count),nk(count),nc(count);
    std::vector<Eigen::Index> input_offset(count+1,0),output_offset(count+1,0);
    for (int item=0;item<count;++item) {
      if (f_values[item]==nullptr || g[item]==nullptr || hg[item]==nullptr ||
          static_cast<int>(f_values[item]->size())!=rank() ||
          static_cast<int>(g[item]->size())!=rank() ||
          static_cast<int>(hg[item]->size())!=rank())
        throw std::invalid_argument(
          "DLR batched prepared convolution: coefficient count != rank");
      nr[item]=(*f_values[item])[0].rows();
      nk[item]=(*f_values[item])[0].cols();
      nc[item]=(*g[item])[0].cols();
      if ((*g[item])[0].rows()!=nk[item] ||
          (*hg[item])[0].rows()!=nk[item] ||
          (*hg[item])[0].cols()!=nc[item])
        throw std::invalid_argument(
          "DLR batched prepared convolution: incompatible matrices");
      for (int l=0;l<rank();++l)
        if ((*f_values[item])[l].rows()!=nr[item] ||
            (*f_values[item])[l].cols()!=nk[item] ||
            (*g[item])[l].rows()!=nk[item] ||
            (*g[item])[l].cols()!=nc[item] ||
            (*hg[item])[l].rows()!=nk[item] ||
            (*hg[item])[l].cols()!=nc[item])
          throw std::invalid_argument(
            "DLR batched prepared convolution: matrix dimensions vary");
      input_offset[item+1]=input_offset[item]+nr[item]*nk[item];
      output_offset[item+1]=output_offset[item]+nr[item]*nc[item];
    }

    MatrixXd f_value_flat(rank(),input_offset.back());
    for (int item=0;item<count;++item)
      for (int l=0;l<rank();++l)
        for (int i=0;i<nr[item];++i)
          for (int j=0;j<nk[item];++j)
            f_value_flat(l,input_offset[item]+i*nk[item]+j)=
              (*f_values[item])[l](i,j);
    const MatrixXd f_flat=values_to_coefficients(f_value_flat);
    MatrixXd hf_flat(rank(),input_offset.back());
    hf_flat.noalias()=thilb_*f_flat;

    MatrixXd diagonal_flat(rank(),output_offset.back());
    MatrixXd off_diagonal_flat(rank(),output_offset.back());
    for (int item=0;item<count;++item) {
      MatrixXd f_l(nr[item],nk[item]),hf_l(nr[item],nk[item]);
      for (int l=0;l<rank();++l) {
        for (int i=0;i<nr[item];++i)
          for (int j=0;j<nk[item];++j) {
            const Eigen::Index column=input_offset[item]+i*nk[item]+j;
            f_l(i,j)=f_flat(l,column);
            hf_l(i,j)=hf_flat(l,column);
          }
        const MatrixXd diagonal=f_l*(*g[item])[l];
        const MatrixXd off_diagonal=
          hf_l*(*g[item])[l]+f_l*(*hg[item])[l];
        for (int i=0;i<nr[item];++i)
          for (int j=0;j<nc[item];++j) {
            const Eigen::Index column=output_offset[item]+i*nc[item]+j;
            diagonal_flat(l,column)=diagonal(i,j);
            off_diagonal_flat(l,column)=off_diagonal(i,j);
          }
      }
    }

    MatrixXd values_flat(rank(),output_offset.back());
    if (reflect_output)
      values_flat.noalias()=beta_*(
        reflected_ttcf2it_*diagonal_flat+
        reflected_cf2it_*off_diagonal_flat);
    else
      values_flat.noalias()=beta_*(
        ttcf2it_*diagonal_flat+cf2it_*off_diagonal_flat);

    std::vector<std::vector<MatrixXd>> values(count);
    for (int item=0;item<count;++item) {
      values[item].assign(rank(),MatrixXd::Zero(nr[item],nc[item]));
      for (int l=0;l<rank();++l)
        for (int i=0;i<nr[item];++i)
          for (int j=0;j<nc[item];++j)
            values[item][l](i,j)=
              values_flat(l,output_offset[item]+i*nc[item]+j);
    }
    return values;
  }

  // Optional caller-owned workspace. No basis-global mutable buffers: MPI
  // workers and independent evaluators never alias this storage. Matrices
  // are reused for repeated shapes; rank-by-M size is bounded by the caller's
  // existing batch cap. No physics data survive a task via this workspace.
  struct ConvolutionScratch {
    bool profile=false,contiguous_products=false;
    double pack_seconds=0,transform_seconds=0,product_seconds=0,unpack_seconds=0;
    std::vector<int> nr,nk,nc;
    std::vector<Eigen::Index> input_offset,output_offset;
    MatrixXd f_value_flat,f_flat,hf_flat,diagonal_flat,off_diagonal_flat,values_flat;
    MatrixXd f_l,hf_l,diagonal,off_diagonal;
    // Each atomic matrix at one DLR node is a contiguous row-major view.
    // Kept separate from the column-major LU buffers to retain the original
    // interpolation solve and its numerical conditioning.
    using RowMatrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
    RowMatrix f_rows,hf_rows,diagonal_rows,off_diagonal_rows;
  };
  std::vector<std::vector<MatrixXd>>
  time_ordered_convolution_values_right_prepared_batch(
      const std::vector<const std::vector<MatrixXd>*>& f_values,
      const std::vector<const std::vector<MatrixXd>*>& g,
      const std::vector<const std::vector<MatrixXd>*>& hg,
      bool reflect_output=false,ConvolutionScratch* scratch=nullptr,
      const std::vector<const std::vector<MatrixXd>*>* left_hilbert=nullptr,
      bool input_coefficients=false,bool output_coefficients=false) const
  {
    // With left_hilbert, f_values already contains coefficients. This lets
    // callers transform a common prefix BEFORE sparse bath branching. Both
    // coefficient arrays then undergo the same time-independent operator.
    if (output_coefficients && (reflect_output || !coefficient_convolution_available()))
      throw std::invalid_argument("DLR coefficient output unavailable or incompatible with reflection");
    if (!scratch && !left_hilbert && !input_coefficients && !output_coefficients)
      return time_ordered_convolution_values_right_prepared_batch_reference(
      f_values,g,hg,reflect_output);
    const int count=static_cast<int>(f_values.size());
    if (static_cast<int>(g.size())!=count ||
        static_cast<int>(hg.size())!=count)
      throw std::invalid_argument(
        "DLR batched prepared convolution: batch sizes differ");
    if (count==0) return {};
    if (left_hilbert && left_hilbert->size()!=f_values.size())
      throw std::invalid_argument("DLR prepared left: batch sizes differ");

    ConvolutionScratch local;
    auto& buffers=scratch ? *scratch : local;
    auto stamp=std::chrono::steady_clock::time_point{};
    const auto mark=[&]() { if (buffers.profile) stamp=std::chrono::steady_clock::now(); };
    const auto elapsed=[&](double& value) {
      if (buffers.profile) value+=std::chrono::duration<double>(std::chrono::steady_clock::now()-stamp).count();
    };
    mark();
    auto& nr=buffers.nr; auto& nk=buffers.nk; auto& nc=buffers.nc;
    nr.resize(count); nk.resize(count); nc.resize(count);
    auto& input_offset=buffers.input_offset; auto& output_offset=buffers.output_offset;
    input_offset.assign(count+1,0); output_offset.assign(count+1,0);
    for (int item=0;item<count;++item) {
      if (f_values[item]==nullptr || g[item]==nullptr || hg[item]==nullptr ||
          static_cast<int>(f_values[item]->size())!=rank() ||
          static_cast<int>(g[item]->size())!=rank() ||
          static_cast<int>(hg[item]->size())!=rank())
        throw std::invalid_argument(
          "DLR batched prepared convolution: coefficient count != rank");
      nr[item]=(*f_values[item])[0].rows();
      nk[item]=(*f_values[item])[0].cols();
      nc[item]=(*g[item])[0].cols();
      if (left_hilbert && (!(*left_hilbert)[item] ||
          (*left_hilbert)[item]->size()!=static_cast<std::size_t>(rank())))
        throw std::invalid_argument("DLR prepared left: coefficient count differs");
      if ((*g[item])[0].rows()!=nk[item] ||
          (*hg[item])[0].rows()!=nk[item] ||
          (*hg[item])[0].cols()!=nc[item])
        throw std::invalid_argument(
          "DLR batched prepared convolution: incompatible matrices");
      for (int l=0;l<rank();++l)
        if ((*f_values[item])[l].rows()!=nr[item] ||
            (*f_values[item])[l].cols()!=nk[item] ||
            (*g[item])[l].rows()!=nk[item] ||
            (*g[item])[l].cols()!=nc[item] ||
            (*hg[item])[l].rows()!=nk[item] ||
            (*hg[item])[l].cols()!=nc[item])
          throw std::invalid_argument(
            "DLR batched prepared convolution: matrix dimensions vary");
      input_offset[item+1]=input_offset[item]+nr[item]*nk[item];
      output_offset[item+1]=output_offset[item]+nr[item]*nc[item];
    }

    auto& f_value_flat=buffers.f_value_flat;
    f_value_flat.resize(rank(),input_offset.back());
    // Flat matrices are column-major: write one contiguous DLR column at a
    // time instead of making rank strided passes across the entire batch.
    for (int item=0;item<count;++item)
      for (int i=0;i<nr[item];++i)
        for (int j=0;j<nk[item];++j)
          for (int l=0;l<rank();++l)
            f_value_flat(l,input_offset[item]+i*nk[item]+j)=
              (*f_values[item])[l](i,j);
    elapsed(buffers.pack_seconds); mark();
    auto& f_flat=buffers.f_flat;
    auto& hf_flat=buffers.hf_flat;
    hf_flat.resize(rank(),input_offset.back());
    if (left_hilbert) {
      f_flat.swap(f_value_flat);
      for (int item=0;item<count;++item) for (int l=0;l<rank();++l) {
        const auto& h=(*(*left_hilbert)[item])[l];
        if (h.rows()!=nr[item] || h.cols()!=nk[item])
          throw std::invalid_argument("DLR prepared left: matrix dimensions differ");
        for (int i=0;i<nr[item];++i) for (int j=0;j<nk[item];++j)
          hf_flat(l,input_offset[item]+i*nk[item]+j)=h(i,j);
      }
    } else {
      if (input_coefficients) f_flat.swap(f_value_flat);
      else f_flat=lu_.solve(f_value_flat);
      hf_flat.noalias()=thilb_*f_flat;
    }

    elapsed(buffers.transform_seconds); mark();
    auto& diagonal_flat=buffers.diagonal_flat;
    auto& off_diagonal_flat=buffers.off_diagonal_flat;
    if (buffers.contiguous_products) {
      using RowMatrix=ConvolutionScratch::RowMatrix;
      // A single blocked layout conversion replaces the repeated strided
      // gather/scatter loops at every atomic product. These include general
      // rectangular nr x nk times nk x nc matrices, not only square pieces.
      auto& f=buffers.f_rows; auto& hf=buffers.hf_rows;
      auto& d=buffers.diagonal_rows; auto& o=buffers.off_diagonal_rows;
      f=f_flat; hf=hf_flat;
      d.resize(rank(),output_offset.back());o.resize(rank(),output_offset.back());
      for (int item=0;item<count;++item) for (int l=0;l<rank();++l) {
        const auto fi=static_cast<Eigen::Index>(l)*input_offset.back()+input_offset[item];
        const auto oi=static_cast<Eigen::Index>(l)*output_offset.back()+output_offset[item];
        Eigen::Map<const RowMatrix> left(f.data()+fi,nr[item],nk[item]);
        Eigen::Map<const RowMatrix> hilbert(hf.data()+fi,nr[item],nk[item]);
        Eigen::Map<RowMatrix> diagonal(d.data()+oi,nr[item],nc[item]);
        Eigen::Map<RowMatrix> off(o.data()+oi,nr[item],nc[item]);
        diagonal.noalias()=left*(*g[item])[l];
        off.noalias()=hilbert*(*g[item])[l]+left*(*hg[item])[l];
      }
      elapsed(buffers.product_seconds);mark();
      auto& values_flat=buffers.values_flat;
      values_flat.resize(rank(),output_offset.back());
      if (output_coefficients)
        values_flat.noalias()=beta_*(ttcf2cf_*d+o);
      else if (reflect_output)
        values_flat.noalias()=beta_*(reflected_ttcf2it_*d+reflected_cf2it_*o);
      else
        values_flat.noalias()=beta_*(ttcf2it_*d+cf2it_*o);
    } else {
      diagonal_flat.resize(rank(),output_offset.back());
      off_diagonal_flat.resize(rank(),output_offset.back());
      for (int item=0;item<count;++item) {
        auto& f_l=buffers.f_l; auto& hf_l=buffers.hf_l;
        f_l.resize(nr[item],nk[item]); hf_l.resize(nr[item],nk[item]);
        auto& diagonal=buffers.diagonal; auto& off_diagonal=buffers.off_diagonal;
        for (int l=0;l<rank();++l) {
          for (int i=0;i<nr[item];++i)
            for (int j=0;j<nk[item];++j) {
              const Eigen::Index column=input_offset[item]+i*nk[item]+j;
              f_l(i,j)=f_flat(l,column);
              hf_l(i,j)=hf_flat(l,column);
            }
          diagonal.noalias()=f_l*(*g[item])[l];
          off_diagonal.noalias()=hf_l*(*g[item])[l]+f_l*(*hg[item])[l];
          for (int i=0;i<nr[item];++i)
            for (int j=0;j<nc[item];++j) {
              const Eigen::Index column=output_offset[item]+i*nc[item]+j;
              diagonal_flat(l,column)=diagonal(i,j);
              off_diagonal_flat(l,column)=off_diagonal(i,j);
            }
        }
      }
  
      elapsed(buffers.product_seconds); mark();
      auto& values_flat=buffers.values_flat;
      values_flat.resize(rank(),output_offset.back());
      // With coefficient output, B*M=T gives h=beta*(M*D+O).
      // The off-diagonal term is ALREADY in the DLR basis. Do not apply B
      // and then solve it away at the next convolution. Matrix products in
      // D and O retain their noncommuting atomic order. The final convolution
      // still uses the existing node/reflected formula at the boundary.
      if (output_coefficients)
        values_flat.noalias()=beta_*(ttcf2cf_*diagonal_flat+off_diagonal_flat);
      else if (reflect_output)
        values_flat.noalias()=beta_*(
          reflected_ttcf2it_*diagonal_flat+
          reflected_cf2it_*off_diagonal_flat);
      else
        values_flat.noalias()=beta_*(
          ttcf2it_*diagonal_flat+cf2it_*off_diagonal_flat);
    }

    elapsed(buffers.transform_seconds); mark();
    std::vector<std::vector<MatrixXd>> values(count);
    for (int item=0;item<count;++item) {
      // Every element below is assigned. Avoid a redundant zero-fill pass.
      values[item].resize(rank());
      for (auto& node:values[item]) node.resize(nr[item],nc[item]);
      for (int i=0;i<nr[item];++i)
        for (int j=0;j<nc[item];++j)
          for (int l=0;l<rank();++l)
            values[item][l](i,j)=
              buffers.values_flat(l,output_offset[item]+i*nc[item]+j);
    }
    elapsed(buffers.unpack_seconds);
    return values;
  }

  // Convert a matrix-valued function between DLR node values and DLR
  // coefficients.  Flattening matrix entries into columns turns the many
  // identical interpolation solves into one dense multiple-right-hand-side
  // solve.  All matrices must have the same dimensions.
  std::vector<MatrixXd>
  matrix_values_to_coefficients(const std::vector<MatrixXd>& values) const
  {
    return transform_matrix_function(values,true);
  }

  std::vector<MatrixXd>
  matrix_coefficients_to_values(
      const std::vector<MatrixXd>& coefficients) const
  {
    return transform_matrix_function(coefficients,false);
  }

  std::vector<MatrixXd>
  time_ordered_convolution_coefficients(
      const std::vector<MatrixXd>& f,
      const std::vector<MatrixXd>& g) const
  {
    return matrix_values_to_coefficients(
      time_ordered_convolution_values(f,g));
  }

private:
  std::vector<MatrixXd>
  transform_matrix_function(const std::vector<MatrixXd>& input,
                            bool values_to_coefficients_transform) const
  {
    if (static_cast<int>(input.size())!=rank())
      throw std::invalid_argument(
        "DLR matrix-function transform: input count != rank");
    if (input.empty()) return {};
    const int rows=input[0].rows(),cols=input[0].cols();
    for (const auto& matrix:input)
      if (matrix.rows()!=rows || matrix.cols()!=cols)
        throw std::invalid_argument(
          "DLR matrix-function transform: matrix dimensions vary");

    MatrixXd flat(rank(),rows*cols);
    for (int l=0;l<rank();++l)
      for (int i=0;i<rows;++i)
        for (int j=0;j<cols;++j)
          flat(l,i*cols+j)=input[l](i,j);
    const MatrixXd transformed=values_to_coefficients_transform
      ? values_to_coefficients(flat) : coefficients_to_values(flat);
    std::vector<MatrixXd> output(
      rank(),MatrixXd::Zero(rows,cols));
    for (int l=0;l<rank();++l)
      for (int i=0;i<rows;++i)
        for (int j=0;j<cols;++j)
          output[l](i,j)=transformed(l,i*cols+j);
    return output;
  }

  void check_rows(Eigen::Index n) const
  {
    if (n!=rank()) throw std::invalid_argument("DLR vector/matrix first dimension != rank");
  }

  void build(int panel_order)
  {
    const auto fine=detail::make_fine_grid(lambda_,panel_order);
    MatrixXd weighted_kernel(fine.time_relative.size(),fine.omega.size());
    for (int i=0; i<weighted_kernel.rows(); ++i)
      for (int j=0; j<weighted_kernel.cols(); ++j)
        weighted_kernel(i,j)=fine.sqrt_time_weights(i)*kernel_relative(fine.time_relative(i),fine.omega(j));

    // Interpolative decomposition of K(t,w): pivot rows of K^T select the
    // frequency skeleton.  These frequencies span the kernel space; they are
    // not fitted poles of any particular G and generally have no direct
    // spectral interpretation.
    auto fpiv=detail::select_rows(weighted_kernel.transpose(),eps_);
    if (fpiv.empty()) throw std::runtime_error("DLR construction returned rank zero");
    std::sort(fpiv.begin(),fpiv.end());
    omega_dimensionless_.resize(fpiv.size());
    for (int i=0; i<static_cast<int>(fpiv.size()); ++i) omega_dimensionless_(i)=fine.omega(fpiv[i]);
    omega_physical_=omega_dimensionless_/beta_;

    MatrixXd time_kernel(fine.time_relative.size(),rank());
    for (int i=0; i<time_kernel.rows(); ++i)
      for (int l=0; l<rank(); ++l)
        time_kernel(i,l)=kernel_relative(fine.time_relative(i),omega_dimensionless_(l));
    // With frequencies fixed, a second pivoting step selects exactly r time
    // nodes.  Their interpolation matrix is square and nonsingular (up to
    // the requested numerical rank), which is the source of the seemingly
    // "magic" reduction from hundreds of tau points to only r nodes.
    auto tpiv=detail::select_rows(time_kernel,0.0,rank());
    std::sort(tpiv.begin(),tpiv.end());
    time_relative_.resize(rank());
    time_physical_.resize(rank());
    cf2it_.resize(rank(),rank());
    for (int i=0; i<rank(); ++i) {
      time_relative_(i)=fine.time_relative(tpiv[i]);
      time_physical_(i)=beta_*relative_to_absolute(time_relative_(i));
      cf2it_.row(i)=time_kernel.row(tpiv[i]);
    }
    lu_.compute(cf2it_);
    if (lu_.matrixLU().diagonal().cwiseAbs().minCoeff()==0.0)
      throw std::runtime_error("DLR interpolation matrix is singular");
    zlu_.compute(cf2it_.cast<std::complex<double> >());

    // R maps values at tau_k directly to values at beta-tau_k:
    //
    //   coefficients = K_nodes^{-1} values,
    //   reflected    = K_reflected coefficients = R values.
    //
    // Forming R once is both faster and exposes reflection as a row-local
    // operation for distributed pair matrices.
    MatrixXd reflected_kernel(rank(),rank());
    for (int i=0; i<rank(); ++i)
      for (int l=0; l<rank(); ++l)
        reflected_kernel(i,l)=kernel_abs(
          1.0-relative_to_absolute(time_relative_(i)),
          omega_dimensionless_(l));
    // Solve R*A=K_ref as A^T*R^T=K_ref^T.  This avoids explicitly forming
    // A^{-1}, whose roundoff is noticeably amplified for a high-accuracy DLR
    // interpolation matrix.
    reflection_values_=cf2it_.transpose().partialPivLu()
      .solve(reflected_kernel.transpose()).transpose();

    // Precompute the divided-difference matrices used by the analytic
    // convolution identity.  thilb_ contains unequal-frequency denominators;
    // ttcf2it_ is the equal-frequency limit proportional to tau*K(tau,w).
    thilb_=MatrixXd::Zero(rank(),rank());
    for (int j=0; j<rank(); ++j)
      for (int k=0; k<rank(); ++k)
        if (j!=k)
          thilb_(j,k)=kernel_abs(0.0,omega_dimensionless_(k))/
                       (omega_dimensionless_(k)-omega_dimensionless_(j));
    ttcf2it_.resize(rank(),rank());
    for (int i=0; i<rank(); ++i) {
      const double tabs=relative_to_absolute(time_relative_(i));
      for (int l=0; l<rank(); ++l)
        ttcf2it_(i,l)=tabs*cf2it_(i,l)*kernel_abs(0.0,omega_dimensionless_(l));
    }
    // A physical diagram needs the second backbone at beta-tau.  Fusing the
    // reflection into its final convolution replaces R*(T*d+K*o) by the two
    // precomputed transforms (R*T)*d and (R*K)*o, avoiding an additional
    // dense DLR transform and a pack/unpack pass for every bath/pole pair.
    reflected_ttcf2it_.resize(rank(),rank());
    reflected_cf2it_.resize(rank(),rank());
    reflected_ttcf2it_.noalias()=reflection_values_*ttcf2it_;
    reflected_cf2it_.noalias()=reflection_values_*cf2it_;
    // One rank-by-rank solve per basis, not an explicit inverse. This adds
    // O(rank^2) storage and no bath- or pseudo-dependent cache. Ordinary
    // convolution formulas above are unchanged when the new flags are off.
    coefficient_convolution_ok_=false;
    // Complete SC bath sums are more sensitive than isolated fixed-bath
    // chains. The validated production regime is 1e-7 or looser; tighter
    // bases retain the node projection at every stage automatically.
    if (eps_>=1.e-7) {
      ttcf2cf_=lu_.solve(ttcf2it_);
      const double residual=(cf2it_*ttcf2cf_-ttcf2it_).norm();
      coefficient_convolution_ok_=ttcf2cf_.allFinite() && std::isfinite(residual) &&
        residual<=1.e-12*std::max(1.e-300,ttcf2it_.norm());
    } else ttcf2cf_.resize(0,0);
  }

  double beta_,omega_max_,eps_,lambda_;
  VectorXd omega_dimensionless_,omega_physical_;
  VectorXd time_relative_,time_physical_;
  MatrixXd cf2it_,reflection_values_,thilb_,ttcf2it_;
  MatrixXd reflected_cf2it_,reflected_ttcf2it_;
  MatrixXd ttcf2cf_;
  bool coefficient_convolution_ok_=false;
  Eigen::PartialPivLU<MatrixXd> lu_;
  Eigen::PartialPivLU<MatrixXcd> zlu_;
};

} // namespace edmft::dlr

#endif
