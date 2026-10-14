#ifndef _VEXP_PROJECT_H_
#define _VEXP_PROJECT_H_

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "Eigen/Dense"

#include <gsl/gsl_multimin.h>

// SharedExpFit implements a variable-projection fit of an N x N matrix-valued
// function Sigma(tau) to a sum of left- and right-decaying exponentials with a
// pole set shared by the whole matrix block:
//
//   Sigma_ab(tau) ~= sum_m Am[m](a,b) * exp(-bm[m] * tau)
//                  + sum_m Bm[m](a,b) * exp(-dm[m] * (beta - tau))
//
// Only the residues Am/Bm depend on the matrix indices (a,b); the exponents
// bm/dm are global for the block. This shared-exponent structure is the key
// property needed later if one wants to keep matrix multiplications in the
// Dyson/Laplace propagation.

class SharedExpFit {
public:
  struct Options {
    int Mleft = 4;
    int Mright = 4;
    int max_iter = 300;
    double tol = 1e-8;
    bool enforce_order = true;
    bool verbose = false;
    Eigen::VectorXd weights;      // optional Nt-sized point weights
    Eigen::VectorXd init_bm;      // optional positive initial guesses
    Eigen::VectorXd init_dm;      // optional positive initial guesses
    double initial_step = 0.5;    // simplex step in log-parameter space
    double min_spacing = 1e-3;    // soft preferred minimum spacing between neighboring poles
    double spacing_penalty = 1e2; // strength of the soft pole-repulsion penalty
  };

  struct Result {
    Eigen::VectorXd bm;                  // shared left exponents   exp(-bm * tau)
    Eigen::VectorXd dm;                  // shared right exponents  exp(-dm * (beta-tau))
    std::vector<Eigen::MatrixXd> Am;     // residues for left exponentials
    std::vector<Eigen::MatrixXd> Bm;     // residues for right exponentials
    double chi2 = std::numeric_limits<double>::infinity();
    int status = GSL_CONTINUE;
    int iterations = 0;
  };

  // tau       : original nonuniform imaginary-time mesh on [0,beta]
  // Sigma_tau : Sigma(tau_i) as one N x N matrix per mesh point
  // beta      : inverse temperature / interval end point
  SharedExpFit(const Eigen::VectorXd& tau,
               const std::vector<Eigen::MatrixXd>& Sigma_tau,
               double beta)
      : SharedExpFit(tau, Sigma_tau, beta, Options()) {}

  SharedExpFit(const Eigen::VectorXd& tau,
               const std::vector<Eigen::MatrixXd>& Sigma_tau,
               double beta,
               const Options& opts)
      : tau_(tau), Sigma_tau_(Sigma_tau), beta_(beta), opts_(opts) {
    if (tau_.size() == 0) throw std::invalid_argument("SharedExpFit: empty tau mesh");
    if (Sigma_tau_.empty()) throw std::invalid_argument("SharedExpFit: empty Sigma_tau");
    if (static_cast<int>(Sigma_tau_.size()) != tau_.size())
      throw std::invalid_argument("SharedExpFit: tau.size() must equal Sigma_tau.size()");
    N_ = static_cast<int>(Sigma_tau_[0].rows());
    if (N_ == 0 || Sigma_tau_[0].cols() != N_)
      throw std::invalid_argument("SharedExpFit: Sigma_tau blocks must be square and non-empty");
    for (const auto& m : Sigma_tau_) {
      if (m.rows() != N_ || m.cols() != N_)
        throw std::invalid_argument("SharedExpFit: all Sigma_tau blocks must have the same shape");
    }
    if (opts_.Mleft < 0 || opts_.Mright < 0)
      throw std::invalid_argument("SharedExpFit: number of exponentials must be non-negative");
    if (opts_.Mleft + opts_.Mright == 0)
      throw std::invalid_argument("SharedExpFit: need at least one exponential");
    if (opts_.weights.size() != 0 && opts_.weights.size() != tau_.size())
      throw std::invalid_argument("SharedExpFit: weights must be empty or tau-sized");
  }

  Result fit() const {
    const int npar = opts_.Mleft + opts_.Mright;
    // GSL sees only the nonlinear parameters (the shared exponents). The
    // linear amplitudes are eliminated analytically inside objective().
    gsl_multimin_function minex_func;
    minex_func.n = npar;
    minex_func.f = &SharedExpFit::objective_gsl;
    minex_func.params = const_cast<SharedExpFit*>(this);

    gsl_vector* x = gsl_vector_alloc(npar);
    gsl_vector* ss = gsl_vector_alloc(npar);
    fill_initial_parameters(x);
    gsl_vector_set_all(ss, opts_.initial_step);

    const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc(T, npar);
    gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

    int status = GSL_CONTINUE;
    int iter = 0;
    while (status == GSL_CONTINUE && iter < opts_.max_iter) {
      ++iter;
      status = gsl_multimin_fminimizer_iterate(s);
      if (status) break;
      const double size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, opts_.tol);
    }

    // Once the minimizer converges, reconstruct the physical quantities
    // (exponents + residue matrices) from the final nonlinear point.
    Result r = evaluate_result(s->x);
    r.status = status;
    r.iterations = iter;

    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(ss);
    gsl_vector_free(x);
    return r;
  }

private:
  struct LinearSolve {
    Eigen::MatrixXd Phi;   // basis matrix evaluated on the tau-mesh
    Eigen::MatrixXd Y;     // flattened matrix-valued input data
    Eigen::MatrixXd Coeff; // linear coefficients for all matrix elements
    double chi2 = std::numeric_limits<double>::infinity();
  };

  Eigen::VectorXd tau_;
  std::vector<Eigen::MatrixXd> Sigma_tau_;
  double beta_;
  Options opts_;
  int N_ = 0;

  static double objective_gsl(const gsl_vector* x, void* params) {
    const auto* self = static_cast<const SharedExpFit*>(params);
    return self->objective(x);
  }

  double objective(const gsl_vector* x) const {
    // Variable projection:
    //   1) convert unconstrained optimizer variables -> positive exponents
    //   2) for these exponents solve the best linear least-squares problem
    //   3) return only the residual norm to the nonlinear optimizer
    const auto [bm, dm] = decode_exponents(x);
    return solve_linear(bm, dm).chi2 + spacing_penalty(bm) + spacing_penalty(dm);
  }

  Result evaluate_result(const gsl_vector* x) const {
    Result r;
    std::tie(r.bm, r.dm) = decode_exponents(x);
    LinearSolve ls = solve_linear(r.bm, r.dm);
    // The linear solve returns one coefficient vector per matrix element.
    // Here we regroup those vectors into one residue matrix per exponent.
    unpack_coefficients(ls.Coeff, r.Am, r.Bm);
    r.chi2 = ls.chi2;
    return r;
  }

  std::pair<Eigen::VectorXd, Eigen::VectorXd> decode_exponents(const gsl_vector* x) const {
    Eigen::VectorXd bm(opts_.Mleft), dm(opts_.Mright);
    int k = 0;
    if (opts_.enforce_order) {
      // Enforce positivity and ordering by accumulating positive increments:
      //   b0 = exp(x0), b1 = b0 + exp(x1), ...
      double accum = 0.0;
      for (int i = 0; i < opts_.Mleft; ++i) {
        accum += std::exp(gsl_vector_get(x, k++));
        bm(i) = accum;
      }
      accum = 0.0;
      for (int i = 0; i < opts_.Mright; ++i) {
        accum += std::exp(gsl_vector_get(x, k++));
        dm(i) = accum;
      }
    } else {
      for (int i = 0; i < opts_.Mleft; ++i) bm(i) = std::exp(gsl_vector_get(x, k++));
      for (int i = 0; i < opts_.Mright; ++i) dm(i) = std::exp(gsl_vector_get(x, k++));
    }
    return {bm, dm};
  }

  void fill_initial_parameters(gsl_vector* x) const {
    Eigen::VectorXd bm0 = opts_.init_bm.size() ? opts_.init_bm : default_left_guess();
    Eigen::VectorXd dm0 = opts_.init_dm.size() ? opts_.init_dm : default_right_guess();
    if (bm0.size() != opts_.Mleft || dm0.size() != opts_.Mright)
      throw std::invalid_argument("SharedExpFit: initial guess size mismatch");

    int k = 0;
    if (opts_.enforce_order) {
      // Convert ordered positive guesses into positive increments expected
      // by decode_exponents().
      double prev = 0.0;
      for (int i = 0; i < opts_.Mleft; ++i) {
        const double delta = (i == 0 ? bm0(i) : bm0(i) - prev);
        gsl_vector_set(x, k++, std::log(std::max(delta, 1e-12)));
        prev = bm0(i);
      }
      prev = 0.0;
      for (int i = 0; i < opts_.Mright; ++i) {
        const double delta = (i == 0 ? dm0(i) : dm0(i) - prev);
        gsl_vector_set(x, k++, std::log(std::max(delta, 1e-12)));
        prev = dm0(i);
      }
    } else {
      for (int i = 0; i < opts_.Mleft; ++i) gsl_vector_set(x, k++, std::log(std::max(bm0(i), 1e-12)));
      for (int i = 0; i < opts_.Mright; ++i) gsl_vector_set(x, k++, std::log(std::max(dm0(i), 1e-12)));
    }
  }

  Eigen::VectorXd default_left_guess() const {
    Eigen::VectorXd out(opts_.Mleft);
    if (opts_.Mleft == 0) return out;
    // Use a simple logarithmic ladder of decay rates as a robust generic
    // starting point. The optimizer can then move these poles where needed.
    const double bmin = 1.0 / std::max(beta_, 1e-12);
    const double bmax = 20.0 / std::max(beta_, 1e-12);
    if (opts_.Mleft == 1) {
      out(0) = std::sqrt(bmin * bmax);
      return out;
    }
    for (int i = 0; i < opts_.Mleft; ++i) {
      const double t = static_cast<double>(i) / static_cast<double>(opts_.Mleft - 1);
      out(i) = std::exp((1.0 - t) * std::log(bmin) + t * std::log(bmax));
    }
    return out;
  }

  Eigen::VectorXd default_right_guess() const {
    Eigen::VectorXd out(opts_.Mright);
    if (opts_.Mright == 0) return out;
    // Right-decaying exponentials are initialized in the same way, but they
    // act on the reflected variable (beta - tau).
    const double dmin = 1.0 / std::max(beta_, 1e-12);
    const double dmax = 20.0 / std::max(beta_, 1e-12);
    if (opts_.Mright == 1) {
      out(0) = std::sqrt(dmin * dmax);
      return out;
    }
    for (int i = 0; i < opts_.Mright; ++i) {
      const double t = static_cast<double>(i) / static_cast<double>(opts_.Mright - 1);
      out(i) = std::exp((1.0 - t) * std::log(dmin) + t * std::log(dmax));
    }
    return out;
  }

  double spacing_penalty(const Eigen::VectorXd& poles) const {
    if (poles.size() < 2 || opts_.spacing_penalty <= 0.0 || opts_.min_spacing <= 0.0)
      return 0.0;

    // Softly penalize neighboring poles that become too close. The optimizer is
    // still free to place nearby poles when strongly supported by the data, but
    // it must "pay" for redundant clustering instead of collapsing many poles
    // into the same narrow region.
    double penalty = 0.0;
    for (int i = 0; i + 1 < poles.size(); ++i) {
      const double gap = poles(i + 1) - poles(i);
      if (gap < opts_.min_spacing) {
        const double deficit = opts_.min_spacing - gap;
        penalty += opts_.spacing_penalty * deficit * deficit;
      }
    }
    return penalty;
  }

  LinearSolve solve_linear(const Eigen::VectorXd& bm, const Eigen::VectorXd& dm) const {
    LinearSolve ls;
    build_basis(bm, dm, ls.Phi);
    build_rhs(ls.Y);

    if (opts_.weights.size()) {
      // Weighted least squares is implemented by premultiplying both the basis
      // and the data rows by sqrt(weight).
      for (int it = 0; it < ls.Phi.rows(); ++it) {
        const double w = std::sqrt(std::max(opts_.weights(it), 0.0));
        ls.Phi.row(it) *= w;
        ls.Y.row(it) *= w;
      }
    }

    // All matrix elements share the same exponents, so we solve one linear
    // least-squares problem with many right-hand sides instead of fitting each
    // (a,b) channel separately. This is exactly the "variable projection"
    // idea: once the poles are fixed, the best residues follow from a purely
    // linear solve.
    ls.Coeff = ls.Phi.colPivHouseholderQr().solve(ls.Y);
    const Eigen::MatrixXd R = ls.Y - ls.Phi * ls.Coeff;
    ls.chi2 = R.squaredNorm();
    return ls;
  }

  void build_basis(const Eigen::VectorXd& bm, const Eigen::VectorXd& dm, Eigen::MatrixXd& Phi) const {
    const int M = opts_.Mleft + opts_.Mright;
    Phi.resize(tau_.size(), M);
    int col = 0;
    // First the exponentials decaying away from tau = 0.
    for (int m = 0; m < opts_.Mleft; ++m, ++col) {
      for (int it = 0; it < tau_.size(); ++it) Phi(it, col) = std::exp(-bm(m) * tau_(it));
    }
    // Then the exponentials decaying away from tau = beta.
    for (int m = 0; m < opts_.Mright; ++m, ++col) {
      for (int it = 0; it < tau_.size(); ++it) Phi(it, col) = std::exp(-dm(m) * (beta_ - tau_(it)));
    }
  }

  void build_rhs(Eigen::MatrixXd& Y) const {
    // Flatten each N x N matrix at fixed tau into one row of a tall matrix:
    //   Y(it, :) = vec( Sigma_tau_[it] )
    // so that Phi * Coeff reconstructs all channels at once.
    const int nch = N_ * N_;
    Y.resize(tau_.size(), nch);
    for (int it = 0; it < tau_.size(); ++it) {
      int col = 0;
      for (int a = 0; a < N_; ++a) {
        for (int b = 0; b < N_; ++b) {
          Y(it, col++) = Sigma_tau_[it](a, b);
        }
      }
    }
  }

  void unpack_coefficients(const Eigen::MatrixXd& C,
                           std::vector<Eigen::MatrixXd>& Am,
                           std::vector<Eigen::MatrixXd>& Bm) const {
    // Undo the vec(N x N) layout and recover one residue matrix per exponent.
    // Row 0..Mleft-1 correspond to left exponentials, the remaining rows to
    // right exponentials.
    const int nch = N_ * N_;
    Am.assign(opts_.Mleft, Eigen::MatrixXd::Zero(N_, N_));
    Bm.assign(opts_.Mright, Eigen::MatrixXd::Zero(N_, N_));
    for (int row = 0; row < C.rows(); ++row) {
      Eigen::MatrixXd M(N_, N_);
      int col = 0;
      for (int a = 0; a < N_; ++a) {
        for (int b = 0; b < N_; ++b) {
          M(a, b) = (col < nch ? C(row, col) : 0.0);
          ++col;
        }
      }
      if (row < opts_.Mleft) Am[row] = M;
      else Bm[row - opts_.Mleft] = M;
    }
  }
};

#endif
