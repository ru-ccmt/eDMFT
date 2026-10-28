// Matrix-valued, symmetry-preserving AAA rational approximation for
// fermionic Matsubara data.  The support-point construction follows
// Appendix A of SciPost Phys. 19, 121 (2025).
#ifndef EDMFT_AAA_H
#define EDMFT_AAA_H

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"
#include "LBFGSpp/include/LBFGSB.h"
#include <iostream>
#include "bblas.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include "aaa_memory.h"

extern "C" {
  void FNAME(zggev)(const char* jobvl, const char* jobvr, const int* n,
                    std::complex<double>* a, const int* lda,
                    std::complex<double>* b, const int* ldb,
                    std::complex<double>* alpha,
                    std::complex<double>* beta,
                    std::complex<double>* vl, const int* ldvl,
                    std::complex<double>* vr, const int* ldvr,
                    std::complex<double>* work, const int* lwork,
                    double* rwork, int* info);
}

namespace edmft::aaa {

using Complex = std::complex<double>;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::VectorXd;

struct Options {
  // Stopping threshold for the maximum Frobenius error over sampled positive
  // Matsubara frequencies.  If relative_tolerance is true it is multiplied
  // by the largest sampled matrix norm.
  double tolerance = 1e-12;
  // One iteration adds a conjugate pair of support points.  The final number
  // of real poles may differ because support points and poles are not the
  // same objects: poles are zeros of the barycentric denominator.
  int max_support_pairs = 25;
  bool relative_tolerance = true;
  // Optional initialization diagnostics. The MPI adapter supplies this only
  // on its master rank; standalone numerical users remain silent by default.
  std::ostream* diagnostics = nullptr;
};

struct PoleExpansion {
  // Final real-axis partial-fraction representation
  //       Delta(z) = sum_l residues[l] / (z-poles[l]).
  // This, rather than the barycentric support representation, is consumed by
  // the diagram evaluator.
  VectorXd poles;
  std::vector<MatrixXcd> residues;

  MatrixXcd evaluate(Complex z) const
  {
    if (residues.empty()) return MatrixXcd();
    MatrixXcd out=MatrixXcd::Zero(residues[0].rows(),residues[0].cols());
    for (int l=0; l<poles.size(); ++l) out.noalias() += residues[l]/(z-poles(l));
    return out;
  }
};

struct Result {
  // z_j and f(z_j) in the AAA barycentric interpolant.  Here z_j are selected
  // Matsubara samples and therefore imaginary; they must not be confused
  // with the (ideally real) zeros of the denominator stored in poles.
  std::vector<Complex> support_points;
  std::vector<MatrixXcd> support_values;
  VectorXcd weights;
  VectorXcd poles;
  std::vector<double> max_error_history;
  double sample_scale = 1.0;

  MatrixXcd evaluate(Complex z) const
  {
    if (support_values.empty()) return MatrixXcd();
    for (int j=0; j<static_cast<int>(support_points.size()); ++j) {
      const double scale=std::max({1.0,std::abs(z),std::abs(support_points[j])});
      if (std::abs(z-support_points[j])<=16.0*std::numeric_limits<double>::epsilon()*scale)
        return support_values[j];
    }
    Complex den=0.0;
    MatrixXcd num=MatrixXcd::Zero(support_values[0].rows(),support_values[0].cols());
    for (int j=0; j<weights.size(); ++j) {
      const Complex q=weights(j)/(z-support_points[j]);
      den += q;
      num.noalias() += q*support_values[j];
    }
    if (std::abs(den)<=std::numeric_limits<double>::min())
      throw std::runtime_error("AAA barycentric denominator vanished at evaluation point");
    return num/den;
  }

  VectorXd real_poles(double imag_tolerance=1e-8) const
  {
    std::vector<double> keep;
    for (int i=0; i<poles.size(); ++i) {
      const double scale=std::max(1.0,std::abs(poles(i).real()));
      if (std::abs(poles(i).imag())<=imag_tolerance*scale) keep.push_back(poles(i).real());
    }
    std::sort(keep.begin(),keep.end());
    VectorXd out(keep.size());
    for (int i=0; i<static_cast<int>(keep.size()); ++i) out(i)=keep[i];
    return out;
  }
};

namespace detail {

inline void validate_samples(const std::vector<Complex>& z,
                             const std::vector<MatrixXcd>& f)
{
  if (z.empty() || z.size()!=f.size()) throw std::invalid_argument("AAA requires equally sized nonempty z and f arrays");
  const int nr=f[0].rows(), nc=f[0].cols();
  if (nr<=0 || nc<=0) throw std::invalid_argument("AAA matrix samples must be nonempty");
  for (int i=0; i<static_cast<int>(z.size()); ++i) {
    if (!(z[i].imag()>0.0) || std::abs(z[i].real())>1e-12*std::max(1.0,std::abs(z[i].imag())))
      throw std::invalid_argument("AAA expects positive, purely imaginary Matsubara points");
    if (f[i].rows()!=nr || f[i].cols()!=nc) throw std::invalid_argument("AAA sample matrices have inconsistent shapes");
  }
}

inline VectorXcd generalized_poles(const std::vector<Complex>& z,
                                   const VectorXcd& w)
{
  const int m=static_cast<int>(z.size());
  if (m<2 || w.size()!=m) return VectorXcd();
  const int n=m+1, lda=n, ldb=n;
  // Linearize sum_j w_j/(z-z_j)=0 as a generalized eigenproblem.  This
  // avoids forming a high-degree polynomial, whose coefficients and roots
  // would be badly conditioned.  Infinite generalized eigenvalues are
  // discarded below.
  MatrixXcd a=MatrixXcd::Zero(n,n), b=MatrixXcd::Zero(n,n);
  for (int j=0; j<m; ++j) {
    a(0,j+1)=w(j);
    a(j+1,0)=1.0;
    a(j+1,j+1)=z[j];
    b(j+1,j+1)=1.0;
  }
  VectorXcd alpha(n), beta(n);
  Complex vl_dummy=0.0, vr_dummy=0.0, work_query=0.0;
  std::vector<double> rwork(8*n);
  const char job='N';
  const int ldv=1;
  int lwork=-1, info=0;
  MatrixXcd aq=a,bq=b;
  FNAME(zggev)(&job,&job,&n,aq.data(),&lda,bq.data(),&ldb,
               alpha.data(),beta.data(),&vl_dummy,&ldv,&vr_dummy,&ldv,
               &work_query,&lwork,rwork.data(),&info);
  if (info!=0) throw std::runtime_error("AAA zggev workspace query failed");
  lwork=std::max(2*n,static_cast<int>(std::ceil(work_query.real())));
  std::vector<Complex> work(lwork);
  aq=a; bq=b;
  FNAME(zggev)(&job,&job,&n,aq.data(),&lda,bq.data(),&ldb,
               alpha.data(),beta.data(),&vl_dummy,&ldv,&vr_dummy,&ldv,
               work.data(),&lwork,rwork.data(),&info);
  if (info!=0) throw std::runtime_error("AAA generalized pole eigenproblem failed, info="+std::to_string(info));

  std::vector<Complex> finite;
  for (int i=0; i<n; ++i) {
    const double threshold=128.0*std::numeric_limits<double>::epsilon()*
                           std::max(1.0,std::abs(alpha(i)));
    if (std::abs(beta(i))>threshold) {
      const Complex p=alpha(i)/beta(i);
      if (std::isfinite(p.real()) && std::isfinite(p.imag())) finite.push_back(p);
    }
  }
  std::sort(finite.begin(),finite.end(),[](Complex a,Complex b) {
    return a.real()==b.real() ? a.imag()<b.imag() : a.real()<b.real();
  });
  VectorXcd out(finite.size());
  for (int i=0; i<static_cast<int>(finite.size()); ++i) out(i)=finite[i];
  return out;
}

struct SymmetricSamples {
  std::vector<Complex> z;
  std::vector<MatrixXcd> f;
};

inline SymmetricSamples make_symmetric_samples(const std::vector<Complex>& zpos,
                                               const std::vector<MatrixXcd>& fpos)
{
  validate_samples(zpos,fpos);
  SymmetricSamples s;
  const int n=zpos.size();
  s.z.reserve(2*n); s.f.reserve(2*n);
  for (int i=0; i<n; ++i) { s.z.push_back(zpos[i]); s.f.push_back(fpos[i]); }
  // Enforcing Delta(z*)=Delta(z)^dagger in the sample set makes the real-pole
  // constraint meaningful and keeps the matrix fit causal/Hermitian up to
  // numerical accuracy.
  for (int i=0; i<n; ++i) { s.z.push_back(std::conj(zpos[i])); s.f.push_back(fpos[i].adjoint()); }
  return s;
}

inline MatrixXcd flatten_samples(const std::vector<MatrixXcd>& f)
{
  const int n=f.size(), nr=f[0].rows(), nc=f[0].cols();
  MatrixXcd out(n,nr*nc);
  for (int k=0; k<n; ++k)
    for (int i=0; i<nr; ++i)
      for (int j=0; j<nc; ++j) out(k,i*nc+j)=f[k](i,j);
  return out;
}

inline std::vector<MatrixXcd> unflatten_residues(const MatrixXcd& c, int nr, int nc)
{
  std::vector<MatrixXcd> out(c.rows(),MatrixXcd(nr,nc));
  for (int l=0; l<c.rows(); ++l)
    for (int i=0; i<nr; ++i)
      for (int j=0; j<nc; ++j) out[l](i,j)=c(l,i*nc+j);
  return out;
}

// A matrix-valued fit can contain many identically zero elements. Their
// optimal residues are zero, so removing these right-hand sides changes
// neither the least-squares objective nor its pole derivative. In particular,
// this is NOT an amplitude cutoff: even a 1e-16 anomalous channel is retained.
// Determine activity AFTER constructing the Hermitian frequency partners;
// this also handles a nonzero off-diagonal entry and its transposed partner.
struct ActiveSamples {
  MatrixXcd values;
  std::vector<int> elements;
  int rows=0, cols=0;

  explicit ActiveSamples(const std::vector<MatrixXcd>& f)
    : rows(f[0].rows()), cols(f[0].cols())
  {
    for (int i=0; i<rows; ++i)
      for (int j=0; j<cols; ++j) {
        bool active=false;
        for (const auto& sample:f)
          if (sample(i,j)!=Complex(0,0)) { active=true; break; }
        if (active) elements.push_back(i*cols+j);
      }
    // Keep one zero RHS for an all-zero input, avoiding special empty-matrix
    // paths in Eigen/L-BFGS-B. This does not introduce a nonzero residue.
    if (elements.empty()) elements.push_back(0);
    values.resize(f.size(),elements.size());
    for (int e=0; e<static_cast<int>(elements.size()); ++e)
      for (int n=0; n<static_cast<int>(f.size()); ++n)
        values(n,e)=f[n](elements[e]/cols,elements[e]%cols);
  }

  std::vector<MatrixXcd> residues(const MatrixXcd& c) const
  {
    std::vector<MatrixXcd> out(c.rows(),MatrixXcd::Zero(rows,cols));
    for (int l=0; l<c.rows(); ++l)
      for (int e=0; e<c.cols(); ++e)
        out[l](elements[e]/cols,elements[e]%cols)=c(l,e);
    return out;
  }
};

inline MatrixXcd design_matrix(const std::vector<Complex>& z, const VectorXd& poles)
{
  MatrixXcd a(z.size(),poles.size());
  for (int n=0; n<a.rows(); ++n)
    for (int l=0; l<a.cols(); ++l) a(n,l)=1.0/(z[n]-poles(l));
  return a;
}

// Variable-projection objective used by the real-pole refinement. QR is
// deliberately retained: forming normal equations would square the condition
// number of the Cauchy matrix. At the optimal residues C, R=B-A*C is orthogonal
// to A, eliminating dC/dp from d||R||/dp.
struct PoleObjective {
  const std::vector<Complex>& z;
  const MatrixXcd& b;
  std::ostream* diagnostics=nullptr;
  int evaluations=0;
  using Clock=std::chrono::steady_clock;
  Clock::time_point start=Clock::now(), last_progress=start;

  double operator()(const VectorXd& p, VectorXd& grad)
  {
    const MatrixXcd a=design_matrix(z,p);
    const MatrixXcd c=a.colPivHouseholderQr().solve(b);
    const MatrixXcd residual=b-a*c;
    const double residual_norm=residual.norm();
    grad.resize(p.size());
    VectorXcd da(z.size());
    for (int l=0; l<p.size(); ++l) {
      for (int n=0; n<static_cast<int>(z.size()); ++n)
        da(n)=1.0/((z[n]-p(l))*(z[n]-p(l)));
      // Contract the frequency index first. Eigen's dot conjugates its LEFT
      // operand, giving sum_n conj(R_ne)*da_n. This is the old gradient
      // exactly, without constructing an n_frequency x n_element matrix dr
      // for each pole. Only the floating-point summation order changes.
      Complex contraction=0;
      for (int e=0; e<c.cols(); ++e)
        contraction+=c(l,e)*residual.col(e).dot(da);
      grad(l)=residual_norm>std::numeric_limits<double>::min()
        ? -contraction.real()/residual_norm : 0.0;
    }
    ++evaluations;
    const auto now=Clock::now();
    if (diagnostics && std::chrono::duration<double>(now-last_progress).count()>=30.) {
      // Evaluations include line-search trials; they are not accepted steps.
      (*diagnostics) << "AAA refinement progress: evaluations=" << evaluations
                     << " elapsed_seconds=" << std::chrono::duration<double>(now-start).count()
                     << " trial_objective=" << residual_norm << std::endl;
      last_progress=now;
    }
    return residual_norm;
  }
};

} // namespace detail

inline Result fit(const std::vector<Complex>& z_positive,
                  const std::vector<MatrixXcd>& f_positive,
                  const Options& opts={})
{
  if (!(opts.tolerance>0.0) || opts.max_support_pairs<1)
    throw std::invalid_argument("AAA tolerance and max_support_pairs must be positive");
  const auto samples=detail::make_symmetric_samples(z_positive,f_positive);
  const int nhalf=z_positive.size(), n=2*nhalf;
  const int nr=f_positive[0].rows(), nc=f_positive[0].cols(), nelem=nr*nc;
  const MatrixXcd fflat=detail::flatten_samples(samples.f);
  MatrixXcd fitflat=MatrixXcd::Zero(n,nelem);
  // Remove ONLY identically zero matrix elements from the least-squares
  // rows. The impurity adapter embeds B independent scalar baths in a BxB
  // diagonal matrix, so B*(B-1) elements would otherwise add nothing but
  // zeros to every Loewner column. Do not use a magnitude threshold: small
  // anomalous channels must survive unchanged. Scan the conjugate-extended
  // samples as well, so a general matrix's transposed partner is retained.
  // Keep the original matrix layout in fflat, fitflat and Result; callers,
  // Hermitian conjugation, norms and channel diagnostics are unaffected.
  std::vector<int> active_elements;
  for (int e=0; e<nelem; ++e)
    if ((fflat.col(e).array()!=Complex(0.0,0.0)).any())
      active_elements.push_back(e);
  // A completely zero input still needs a nonempty homogeneous SVD problem
  // to construct the (zero-valued) barycentric result.
  if (active_elements.empty()) active_elements.push_back(0);
  const Eigen::Index nactive=active_elements.size();
  // Reserve one maximum-sized backing buffer for the whole greedy search.
  // Reallocating progressively larger matrices at every support count can
  // leave their freed pages retained by the system allocator until much
  // later. A block view changes only the logical shape, never the allocation.
  const int max_support=2*std::min(opts.max_support_pairs,nhalf);
  MatrixXd real_workspace(2*static_cast<Eigen::Index>(std::max(n-2,0))*nactive,
                          max_support);
  if (opts.diagnostics) {
    *opts.diagnostics << "AAA fit workspace: positive_samples=" << nhalf
      << " matrix_elements=" << nelem << " active_elements=" << nactive
      << " construction=direct-real-reused-inplace-qr"
      << " reserved_mib=" << double(real_workspace.size())*sizeof(double)/(1024.0*1024.0)
      << std::endl;
    detail::print_memory(opts.diagnostics,"fit-start");
  }
  std::vector<bool> chosen(n,false);
  std::vector<int> support, remaining;
  Result result;
  result.sample_scale=0.0;
  for (const auto& f:f_positive) result.sample_scale=std::max(result.sample_scale,f.norm());
  result.sample_scale=std::max(result.sample_scale,1.0);

  for (int pair=0; pair<opts.max_support_pairs; ++pair) {
    // Greedy AAA step: select the as-yet unused positive Matsubara sample
    // with the largest current matrix error and add its conjugate partner at
    // the same time.  Barycentric interpolation will be exact at both.
    int best=-1; double besterr=-1.0;
    for (int i=0; i<nhalf; ++i) if (!chosen[i]) {
      const double err=(fflat.row(i)-fitflat.row(i)).norm();
      if (err>besterr) { besterr=err; best=i; }
    }
    if (best<0) break;
    support.push_back(best); support.push_back(best+nhalf);
    chosen[best]=chosen[best+nhalf]=true;
    remaining.clear();
    for (int i=0; i<n; ++i) if (!chosen[i]) remaining.push_back(i);

    // The Loewner matrix has one row per unused sample and matrix element,
    // but only m columns.  It is tall, not rank-deficient merely because it
    // is rectangular: only min(rows,columns)=m singular values exist.  The
    // right singular vector belonging to the smallest singular value gives
    // the unit-norm weights that minimize ||A w||_2.
    const int m=support.size(), k=m/2;
    const Eigen::Index nrows=static_cast<Eigen::Index>(remaining.size())*nactive;
    // Parameterize conjugate weights by k complex numbers and rewrite the
    // complex least-squares problem as a real one.  This imposes the symmetry
    // exactly rather than hoping an unconstrained SVD rediscovers it.
    auto areal=real_workspace.topLeftCorner(2*nrows,2*k);
    // Form Re/Im of A_plus+A_minus and i*(A_plus-A_minus) directly.
    // The former implementation retained both the complex Loewner matrix
    // and a second complex transformed copy alongside areal and QR storage.
    // Neither is needed. Column-pair-major traversal writes contiguous
    // columns, with exactly the same conjugate-weight parameterization.
    for (int j=0; j<k; ++j) {
      const int sp=support[2*j], sm=support[2*j+1];
      for (Eigen::Index q=0; q<static_cast<Eigen::Index>(remaining.size()); ++q) {
        const int ii=remaining[q];
        const Complex dp=samples.z[ii]-samples.z[sp];
        const Complex dm=samples.z[ii]-samples.z[sm];
        for (Eigen::Index a=0; a<nactive; ++a) {
          const int e=active_elements[a];
          const Complex ap=(fflat(ii,e)-fflat(sp,e))/dp;
          const Complex am=(fflat(ii,e)-fflat(sm,e))/dm;
          const Complex left=ap+am, right=Complex(0.0,1.0)*(ap-am);
          const Eigen::Index row=q*nactive+a;
          areal(row,j)=left.real();
          areal(nrows+row,j)=left.imag();
          areal(row,k+j)=right.real();
          areal(nrows+row,k+j)=right.imag();
        }
      }
    }
    if (opts.diagnostics) {
      *opts.diagnostics << "AAA fit SVD: support_points=" << m
        << " rows=" << areal.rows() << " columns=" << m
        << " matrix_mib=" << double(areal.size())*sizeof(double)/(1024.0*1024.0)
        << std::endl;
      detail::print_memory(opts.diagnostics,"before-svd");
    }
    // Keep the stable pivoted QR -> small SVD solve, but perform the tall QR
    // IN PLACE. Eigen's direct JacobiSVD otherwise owns a second tall copy,
    // whose shape changes at each support count. For A*P=Q*R, the minimizing
    // right vector of A is P*v_min(R); no Q or left singular vectors are
    // needed. Never form A^T*A, which squares the condition number.
    VectorXd wr;
    if (areal.rows()>=m) {
      const double scale=areal.cwiseAbs().maxCoeff();
      if (!std::isfinite(scale)) throw std::runtime_error("nonfinite AAA Loewner matrix");
      if (scale>0.0) areal/=scale;
      Eigen::ColPivHouseholderQR<Eigen::Ref<MatrixXd>> qr(areal);
      const MatrixXd r=qr.matrixQR().topLeftCorner(m,m).template triangularView<Eigen::Upper>();
      Eigen::JacobiSVD<MatrixXd> svd(r,Eigen::ComputeFullV);
      wr=qr.colsPermutation()*svd.matrixV().col(m-1);
    } else if (areal.rows()>0) {
      // Very small data sets can have fewer unused rows than supports.
      // Full V includes the nullspace in this wide-matrix corner case.
      Eigen::JacobiSVD<MatrixXd> svd(areal,Eigen::ComputeFullV);
      wr=svd.matrixV().col(m-1);
    } else {
      // All samples are now supports: any nonzero weights interpolate them.
      wr=VectorXd::Zero(m); wr.head(k).setOnes();
    }
    detail::print_memory(opts.diagnostics,"after-svd");
    result.weights.resize(m);
    for (int j=0; j<k; ++j) {
      result.weights(2*j)=Complex(wr(j),wr(k+j));
      result.weights(2*j+1)=std::conj(result.weights(2*j));
    }

    result.support_points.resize(m);
    result.support_values.resize(m);
    for (int j=0; j<m; ++j) {
      result.support_points[j]=samples.z[support[j]];
      result.support_values[j]=samples.f[support[j]];
      fitflat.row(support[j])=fflat.row(support[j]);
    }
    for (const int ii:remaining) {
      Complex den=0.0;
      Eigen::RowVectorXcd num=Eigen::RowVectorXcd::Zero(nelem);
      for (int j=0; j<m; ++j) {
        const Complex q=result.weights(j)/(samples.z[ii]-result.support_points[j]);
        den+=q; num.noalias()+=q*fflat.row(support[j]);
      }
      fitflat.row(ii)=num/den;
    }
    double err=0.0;
    for (int i=0; i<nhalf; ++i) err=std::max(err,(fitflat.row(i)-fflat.row(i)).norm());
    result.max_error_history.push_back(err);
    const double threshold=opts.tolerance*(opts.relative_tolerance ? result.sample_scale : 1.0);
    if (err<=threshold || static_cast<int>(support.size())>=n) break;
  }
  result.poles=detail::generalized_poles(result.support_points,result.weights);
  detail::print_memory(opts.diagnostics,"fit-done");
  return result;
}

inline PoleExpansion fit_residues(const std::vector<Complex>& z_positive,
                                  const std::vector<MatrixXcd>& f_positive,
                                  const VectorXd& poles)
{
  if (poles.size()==0) return {};
  const auto samples=detail::make_symmetric_samples(z_positive,f_positive);
  const detail::ActiveSamples active(samples.f);
  const MatrixXcd a=detail::design_matrix(samples.z,poles);
  // Once pole locations are fixed the residues enter linearly.  Solve all
  // matrix elements together; every column shares the same Cauchy matrix A.
  const MatrixXcd c=a.colPivHouseholderQr().solve(active.values);
  PoleExpansion out;
  out.poles=poles;
  out.residues=active.residues(c);
  // With real poles and conjugate-paired samples the residues are Hermitian.
  // Enforce that exact structure at roundoff level.
  if (f_positive[0].rows()==f_positive[0].cols())
    for (auto& r:out.residues) r=0.5*(r+r.adjoint()).eval();
  return out;
}

struct RefinementReport {
  PoleExpansion expansion;
  int iterations=0;
  int evaluations=0;
  double seconds=0.0;
  double objective=0.0;
  double relative_rms_error=0.0;
  double max_error=0.0;
};

// Reduce an already accurate real-pole expansion to a prescribed number of
// poles.  This is intended for higher-order diagrams, whose cost grows as the
// square of the pole count.  At every step we try removing each remaining
// pole, refit all residues by linear least squares, and retain the deletion
// with the smallest worst-channel relative RMS error.  The minimax channel
// score is important in a superconducting calculation: a small anomalous
// channel must not be sacrificed merely because the normal channels dominate
// the global Frobenius norm.
struct CompressionReport {
  PoleExpansion expansion;
  std::vector<double> channel_relative_rms;
  double relative_rms_error=0.0;
  double worst_channel_relative_rms=0.0;
};

inline CompressionReport compress_poles_greedy(
    const std::vector<Complex>& z_positive,
    const std::vector<MatrixXcd>& f_positive,
    const PoleExpansion& source,
    int target_poles,
    double channel_scale_floor_ratio=1.e-12,
    bool preserve_both_frequency_signs=false)
{
  detail::validate_samples(z_positive,f_positive);
  if (source.poles.size()==0)
    throw std::invalid_argument("AAA pole compression requires a nonempty source expansion");
  if (target_poles<1 || target_poles>source.poles.size())
    throw std::invalid_argument("AAA compressed pole count is outside [1,source pole count]");

  const int channels=f_positive[0].rows();
  if (f_positive[0].cols()!=channels)
    throw std::invalid_argument("AAA channel-preserving compression expects square samples");
  std::vector<double> scale2(channels,0.0);
  for (const auto& value:f_positive)
    for (int channel=0;channel<channels;++channel)
      scale2[channel]+=std::norm(value(channel,channel));
  const double largest_scale2=*std::max_element(scale2.begin(),scale2.end());
  const double scale2_floor=largest_scale2*channel_scale_floor_ratio*
                            channel_scale_floor_ratio;

  // Production hybridization samples are diagonal in the independent CIX
  // channel representation.  Solve only those physical columns here instead
  // of flattening the auxiliary channels-by-channels matrix: for 12 channels
  // this reduces residue-refit work by another factor of twelve.
  const auto symmetric=detail::make_symmetric_samples(
    z_positive,f_positive);
  MatrixXcd diagonal_samples(
    symmetric.f.size(),channels);
  double off_diagonal_norm2=0.0,total_sample_norm2=0.0;
  for (int n=0;n<diagonal_samples.rows();++n) {
    total_sample_norm2+=symmetric.f[n].squaredNorm();
    for (int channel=0;channel<channels;++channel) {
      diagonal_samples(n,channel)=symmetric.f[n](channel,channel);
      off_diagonal_norm2+=
        symmetric.f[n].row(channel).squaredNorm()-
        std::norm(symmetric.f[n](channel,channel));
    }
  }
  const bool diagonal_input=off_diagonal_norm2<=
    256.0*std::numeric_limits<double>::epsilon()*
    std::max(total_sample_norm2,std::numeric_limits<double>::min());
  const auto fit_candidate=[&](const VectorXd& candidate_poles) {
    if (!diagonal_input)
      return fit_residues(z_positive,f_positive,candidate_poles);
    const MatrixXcd design=detail::design_matrix(
      symmetric.z,candidate_poles);
    const MatrixXcd coefficients=
      design.colPivHouseholderQr().solve(diagonal_samples);
    PoleExpansion expansion;
    expansion.poles=candidate_poles;
    expansion.residues.assign(
      candidate_poles.size(),MatrixXcd::Zero(channels,channels));
    for (int pole=0;pole<candidate_poles.size();++pole)
      for (int channel=0;channel<channels;++channel)
        expansion.residues[pole](channel,channel)=
          coefficients(pole,channel).real();
    return expansion;
  };

  const auto errors=[&](const PoleExpansion& expansion) {
    CompressionReport report;
    report.expansion=expansion;
    report.channel_relative_rms.assign(channels,-1.0);
    std::vector<double> error2(channels,0.0);
    double total_error2=0.0,total_scale2=0.0;
    for (int n=0;n<static_cast<int>(z_positive.size());++n) {
      const MatrixXcd approximation=expansion.evaluate(z_positive[n]);
      const MatrixXcd difference=approximation-f_positive[n];
      total_error2+=difference.squaredNorm();
      total_scale2+=f_positive[n].squaredNorm();
      for (int channel=0;channel<channels;++channel)
        error2[channel]+=std::norm(difference(channel,channel));
    }
    report.relative_rms_error=std::sqrt(
      total_error2/std::max(total_scale2,std::numeric_limits<double>::min()));
    bool found=false;
    for (int channel=0;channel<channels;++channel)
      if (scale2[channel]>std::max(scale2_floor,
                                  std::numeric_limits<double>::min())) {
        report.channel_relative_rms[channel]=
          std::sqrt(error2[channel]/scale2[channel]);
        report.worst_channel_relative_rms=std::max(
          report.worst_channel_relative_rms,
          report.channel_relative_rms[channel]);
        found=true;
      }
    if (!found)
      report.worst_channel_relative_rms=
        std::numeric_limits<double>::infinity();
    return report;
  };

  VectorXd poles=source.poles;
  while (poles.size()>target_poles) {
    bool have_best=false;
    CompressionReport best;
    VectorXd best_poles;
    for (int removed=0;removed<poles.size();++removed) {
      // FCA can request one pole on each side of zero. This optional guard
      // never changes the existing TCA compression policy. It protects both
      // signs only when present in the source and when the budget permits.
      if (preserve_both_frequency_signs && target_poles>=2 &&
          (poles.array()<0.0).any() && (poles.array()>0.0).any()) {
        if (poles(removed)<0.0 && (poles.array()<0.0).count()==1) continue;
        if (poles(removed)>0.0 && (poles.array()>0.0).count()==1) continue;
      }
      VectorXd candidate_poles(poles.size()-1);
      for (int source_index=0,target_index=0;
           source_index<poles.size();++source_index)
        if (source_index!=removed)
          candidate_poles(target_index++)=poles(source_index);
      CompressionReport candidate=errors(fit_candidate(candidate_poles));
      const bool better=!have_best ||
        candidate.worst_channel_relative_rms<
          best.worst_channel_relative_rms ||
        (candidate.worst_channel_relative_rms==
           best.worst_channel_relative_rms &&
         candidate.relative_rms_error<best.relative_rms_error);
      if (better) {
        have_best=true;
        best=std::move(candidate);
        best_poles=std::move(candidate_poles);
      }
    }
    poles=std::move(best_poles);
  }
  return errors(fit_candidate(poles));
}

inline RefinementReport refine_poles(const std::vector<Complex>& z_positive,
                                     const std::vector<MatrixXcd>& f_positive,
                                     VectorXd initial_poles,
                                     double lower_bound,
                                     double upper_bound,
                                     int max_iterations=200,
                                     double gradient_tolerance=1e-10,
                                     std::ostream* diagnostics=nullptr)
{
  if (initial_poles.size()==0) return {};
  if (!(lower_bound<upper_bound)) throw std::invalid_argument("AAA pole bounds must be ordered");
  const auto start=std::chrono::steady_clock::now();
  detail::print_memory(diagnostics,"refinement-start");
  const auto samples=detail::make_symmetric_samples(z_positive,f_positive);
  const detail::ActiveSamples active(samples.f);
  if (diagnostics)
    (*diagnostics) << "AAA refinement start: poles=" << initial_poles.size()
                   << " samples=" << samples.z.size()
                   << " active_elements=" << active.values.cols()
                   << " full_elements=" << active.rows*active.cols << std::endl;

  // AAA itself does not guarantee real poles.  We retain the real candidates
  // and optimize their real positions with box-constrained L-BFGS-B.  For
  // each trial pole vector the optimal residues are eliminated by a linear
  // least-squares solve (variable projection), so the nonlinear optimizer
  // sees only one scalar variable per pole.
  detail::PoleObjective objective{samples.z,active.values,diagnostics};

  LBFGSpp::LBFGSBParam<double> param;
  param.epsilon=gradient_tolerance;
  param.epsilon_rel=gradient_tolerance;
  param.max_iterations=max_iterations;
  LBFGSpp::LBFGSBSolver<double> solver(param);
  VectorXd lb=VectorXd::Constant(initial_poles.size(),lower_bound);
  VectorXd ub=VectorXd::Constant(initial_poles.size(),upper_bound);
  double fx=0.0;
  const int iterations=solver.minimize(objective,initial_poles,fx,lb,ub);
  std::sort(initial_poles.data(),initial_poles.data()+initial_poles.size());

  RefinementReport report;
  report.expansion=fit_residues(z_positive,f_positive,initial_poles);
  report.iterations=iterations;
  report.evaluations=objective.evaluations;
  report.objective=fx;
  double scale2=0.0,err2=0.0,maxerr=0.0;
  for (int i=0; i<static_cast<int>(z_positive.size()); ++i) {
    const MatrixXcd d=report.expansion.evaluate(z_positive[i])-f_positive[i];
    err2+=d.squaredNorm(); scale2+=f_positive[i].squaredNorm();
    maxerr=std::max(maxerr,d.norm());
  }
  report.relative_rms_error=std::sqrt(err2/std::max(scale2,std::numeric_limits<double>::min()));
  report.max_error=maxerr;
  report.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  if (diagnostics)
    (*diagnostics) << "AAA refinement done: iterations=" << report.iterations
                   << " evaluations=" << report.evaluations
                   << " seconds=" << report.seconds
                   << " relative_rms=" << report.relative_rms_error << std::endl;
  detail::print_memory(diagnostics,"refinement-done");
  return report;
}

} // namespace edmft::aaa

#include "aaa_symmetry.h"

#endif
