#ifndef EDMFT_ANALYTIC_ORDERS_H
#define EDMFT_ANALYTIC_ORDERS_H

#include <cmath>
#include <stdexcept>
#include <string>

namespace edmft::orders {

// N_min_order is the analytic normalization anchor, NOT the first order
// absent from the analytic result. MC visits [N_min_order,Norder], using
// its lowest sector only for normalization. InterpolateOnRightMesh adds
// sectors strictly above that anchor (plus the existing cOCA correction).
// Keep this routing shared by workspace initialization and the iteration
// driver: otherwise legacy OCA + analytic TCA + MC can miss its DLR setup.
struct Plan {
  bool tca = false;
  bool fca = false;
  bool sample_higher_orders = false;
};

inline Plan make_plan(bool imaginary_axis, int minimum, int maximum,
                      const std::string& tca_solver,
                      const std::string& fca_solver)
{
  if (minimum<1 || maximum<minimum)
    throw std::invalid_argument("require 1 <= N_min_order <= Norder");
  if (minimum>4 || (!imaginary_axis && minimum>2))
    throw std::invalid_argument(
      "N_min_order has no analytic normalization: supported anchors are "
      "NCA=1, OCA=2, and on the imaginary axis TCA=3, FCA=4");
  Plan result;
  result.sample_higher_orders = maximum>minimum;
  result.tca = imaginary_axis && minimum>=3;
  result.fca = imaginary_axis && minimum>=4;
  if (result.tca && tca_solver!="dlr")
    throw std::invalid_argument(
      "N_min_order>=3 requires TCA_solver=dlr for the analytic reference; "
      "to sample TCA with MC use N_min_order=1 or 2");
  if (result.fca && fca_solver!="dlr")
    throw std::invalid_argument(
      "N_min_order=4 requires FCA_solver=dlr for the analytic reference; "
      "to sample FCA with MC use N_min_order<=3");
  return result;
}

// Phi can be signed: a negative reference is not in itself an error. Zero
// cannot normalize signed MC measurements, and NaN/Inf must never propagate
// into an apparently successful output. Do not introduce an arbitrary
// small-Phi cutoff, since physically small nonzero references are valid.
inline bool valid_reference(double phi)
{
  return std::isfinite(phi) && phi!=0.0;
}

} // namespace edmft::orders
#endif
