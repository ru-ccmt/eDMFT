// Common compact-imaginary-axis data used by OCA and higher-order diagrams.
#ifndef EDMFT_IMAG_AXIS_DLR_H
#define EDMFT_IMAG_AXIS_DLR_H

#include "aaa.h"
#include "dlr.h"
#include "symmetry_rows.h"
#include "oca_symmetry.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace edmft::imag_dlr {

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Numerical data shared by every compact imaginary-axis diagram evaluator.
//
// There are two independent compact representations in this object:
//
//   * AAA supplies real physical pole positions and channel residues for the
//     hybridization.  Its poles are not the DLR frequencies.
//   * DLR supplies a fixed interpolation basis for all pseudo propagators and
//     for all intermediate convolution results.
//
// Keeping both in one workspace is useful because products G(tau)K(tau,w_l)
// occur in every OCA cut and will occur again at third order.  They can be
// formed at only r DLR time nodes and transformed back to r coefficients.
class Workspace {
public:
  Workspace() = default;
  void set_spatial_symmetry(oca_symmetry::Verification map) { spatial_symmetry_=std::move(map); }
  const oca_symmetry::Verification& spatial_symmetry() const { return spatial_symmetry_; }
  void set_self_rows(time_reversal::Rows rows) { self_rows_=std::move(rows); }
  const time_reversal::Rows& self_rows() const {return self_rows_;}
  void set_atomic_optimizations(bool subblocks,bool skip_zero_roots) {
    atomic_subblocks_=subblocks;atomic_skip_zero_roots_=skip_zero_roots;
  }
  bool atomic_subblocks() const {return atomic_subblocks_;}
  bool atomic_skip_zero_roots() const {return atomic_skip_zero_roots_;}

  void prepare_hybridization(const dlr::FermionicBasis& basis,
                             const aaa::PoleExpansion& expansion,
                             int channel_count,
                             double imaginary_tolerance=1.e-9)
  {
    if (expansion.poles.size()==0 || expansion.residues.empty())
      throw std::invalid_argument("imaginary-axis DLR workspace requires a nonempty AAA expansion");
    if (static_cast<int>(expansion.residues.size())!=expansion.poles.size())
      throw std::invalid_argument("AAA pole/residue count mismatch");
    if (channel_count<=0) throw std::invalid_argument("hybridization channel count must be positive");

    // The basis is owned by ImagAxis0 and outlives this workspace.  Storing a
    // pointer avoids copying its interpolation matrices and LU factors.
    self_rows_={}; spatial_symmetry_={};
    basis_=&basis;
    poles_=expansion.poles;
    const int r=basis.rank(), np=poles_.size();
    // The current solver diagonalizes hybridization channels before entering
    // diagram evaluation.  Hence only the diagonal of each matrix residue is
    // needed here.  Rejecting an appreciable imaginary part catches a broken
    // real-axis pole fit rather than silently discarding it.
    residues_.resize(channel_count,np);
    for (int l=0; l<np; ++l) {
      const auto& residue=expansion.residues[l];
      if (residue.rows()<channel_count || residue.cols()<channel_count)
        throw std::invalid_argument("AAA residue matrix is smaller than the hybridization channel space");
      for (int s=0; s<channel_count; ++s) {
        const auto z=residue(s,s);
        if (std::abs(z.imag())>imaginary_tolerance*std::max(1.0,std::abs(z.real())))
          throw std::runtime_error("imaginary-axis OCA requires real channel residues");
        residues_(s,l)=z.real();
      }
    }

    // Precompute both signs because reflection of a line changes w to -w.
    // kernel_abs is evaluated in the overflow-safe branch appropriate to the
    // sign, so these tables stay O(1) even for beta*|omega| >> 1.
    // Rows are DLR time nodes; columns are AAA poles.
    kernel_plus_.resize(r,np);
    kernel_minus_.resize(r,np);
    kernel_zero_plus_.resize(np);
    kernel_zero_minus_.resize(np);
    for (int l=0; l<np; ++l) {
      const double w=beta()*poles_(l);
      kernel_zero_plus_(l)=dlr::kernel_abs(0.0,w);
      kernel_zero_minus_(l)=dlr::kernel_abs(0.0,-w);
      for (int k=0; k<r; ++k) {
        const double t=basis.time_nodes()(k)/beta();
        kernel_plus_(k,l)=dlr::kernel_abs(t,w);
        kernel_minus_(k,l)=dlr::kernel_abs(t,-w);
      }
    }
    // These are standard spectral-kernel sums at the DLR nodes.  Note that
    // the legacy FourierOfDelta object uses the opposite overall sign.  The
    // physical-cut evaluator documents and applies the one remaining sign;
    // in self-energy diagrams two hybridizations make the signs cancel.
    delta_plus_=kernel_plus_*residues_.transpose();
    delta_minus_=kernel_minus_*residues_.transpose();
  }

  void update_pseudos(const MatrixXd& coefficients)
  {
    require_ready();
    if (coefficients.rows()!=rank())
      throw std::invalid_argument("pseudo-particle DLR coefficient count differs from workspace rank");
    // Column p contains the r DLR coefficients of pseudo propagator G_p.
    // Cache node values too: pointwise kernel dressing is most stable in
    // value space and the same values are reused for every AAA pole.
    g_coefficients_=coefficients;
    g_values_=basis_->coefficients_to_values(coefficients);
    // The undressed propagator is also the most common right-hand factor in
    // analytic convolutions.  Cache H*G once so callers using the prepared
    // convolution kernel need no repeated r-by-r transform for this base
    // case.  Pole-dressed variants are prepared locally by their evaluator.
    g_hilbert_coefficients_=basis_->convolution_hilbert(coefficients);
  }

  bool ready() const
  {
    return basis_ && poles_.size()>0 && residues_.size()>0;
  }
  bool has_pseudos() const { return g_coefficients_.size()>0; }
  int rank() const { require_ready(); return basis_->rank(); }
  int pole_count() const { return poles_.size(); }
  int channel_count() const { return residues_.rows(); }
  int pseudo_count() const { return g_coefficients_.cols(); }
  double beta() const { require_ready(); return basis_->beta(); }

  const dlr::FermionicBasis& basis() const { require_ready(); return *basis_; }
  const VectorXd& poles() const { return poles_; }
  const MatrixXd& residues() const { return residues_; }
  const MatrixXd& kernel_plus() const { return kernel_plus_; }
  const MatrixXd& kernel_minus() const { return kernel_minus_; }
  const VectorXd& kernel_zero_plus() const { return kernel_zero_plus_; }
  const VectorXd& kernel_zero_minus() const { return kernel_zero_minus_; }
  const MatrixXd& delta_plus() const { return delta_plus_; }
  const MatrixXd& delta_minus() const { return delta_minus_; }
  const MatrixXd& pseudo_coefficients() const { return g_coefficients_; }
  const MatrixXd& pseudo_hilbert_coefficients() const
  {
    return g_hilbert_coefficients_;
  }
  const MatrixXd& pseudo_values() const { return g_values_; }

  // Coefficients of G_p(tau) K(tau,+/-omega_l), with one output column for
  // every pseudo-particle index in pseudos.
  //
  // We multiply at the r interpolation nodes and solve the already factored
  // interpolation system.  Never split K into exp(-tau*w) and a Fermi factor:
  // either piece can overflow although their product is perfectly bounded.
  // frequency_sign selects the stable K(tau,+w_l) or K(tau,-w_l) table; it
  // does not change the AAA pole or residue stored in the workspace.
  MatrixXd tilted_pseudo_coefficients(const std::vector<int>& pseudos,
                                      int pole,
                                      int frequency_sign) const
  {
    require_pseudos();
    if (pole<0 || pole>=pole_count()) throw std::out_of_range("AAA pole index");
    const auto kernel=(frequency_sign>=0 ? kernel_plus_.col(pole) : kernel_minus_.col(pole));
    MatrixXd values(rank(),pseudos.size());
    for (int j=0; j<static_cast<int>(pseudos.size()); ++j) {
      const int p=pseudos[j];
      if (p<0 || p>=pseudo_count()) throw std::out_of_range("pseudo-particle index");
      values.col(j)=g_values_.col(p).array()*kernel.array();
    }
    return basis_->values_to_coefficients(values);
  }

private:
  time_reversal::Rows self_rows_;
  oca_symmetry::Verification spatial_symmetry_;
  void require_ready() const
  {
    if (!ready()) throw std::logic_error("imaginary-axis DLR workspace is not prepared");
  }
  void require_pseudos() const
  {
    require_ready();
    if (!has_pseudos()) throw std::logic_error("imaginary-axis DLR workspace has no pseudo-particles");
  }

  bool atomic_subblocks_=true,atomic_skip_zero_roots_=true;
  const dlr::FermionicBasis* basis_=nullptr;
  VectorXd poles_, kernel_zero_plus_, kernel_zero_minus_;
  MatrixXd residues_, kernel_plus_, kernel_minus_;
  MatrixXd delta_plus_, delta_minus_;
  MatrixXd g_coefficients_, g_values_, g_hilbert_coefficients_;
};

} // namespace edmft::imag_dlr

#endif
