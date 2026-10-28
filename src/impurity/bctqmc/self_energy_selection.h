#ifndef EDMFT_SELF_ENERGY_SELECTION_H
#define EDMFT_SELF_ENERGY_SELECTION_H
#include <stdexcept>
#include <string>
#include "json.hpp"
namespace edmft::self_energy {
enum class Estimator { dyson, bulla, kugler, blended };
inline Estimator selected(const nlohmann::json& pars) {
  const auto it=pars.find("self_energy_estimator");
  if (it!=pars.end() && !it->is_string())
    throw std::invalid_argument("self_energy_estimator must be dyson, bulla, kugler, or blended");
  const std::string name=it==pars.end() ? "legacy" : it->get<std::string>();
  if (name=="dyson") return Estimator::dyson;
  if (name=="bulla") return Estimator::bulla;
  if (name=="kugler") return Estimator::kugler;
  if (name=="blended") return Estimator::blended;
  // Omitted selector preserves historical boolean and alpha settings.
  if (name=="legacy") return pars.value("Sdyson_Bulla_blend",true)
    ? Estimator::blended : Estimator::bulla;
  throw std::invalid_argument("self_energy_estimator must be dyson, bulla, kugler, or blended");
}
inline const char* name(Estimator e) {
  switch(e) {
  case Estimator::dyson: return "dyson";
  case Estimator::bulla: return "bulla";
  case Estimator::kugler: return "kugler";
  case Estimator::blended: return "blended";
  }
  throw std::logic_error("invalid self-energy estimator");
}
inline void validate(const nlohmann::json& pars,bool imaginary_axis,bool kugler_available) {
  const auto e=selected(pars);
  if (!imaginary_axis && pars.value("self_energy_estimator",std::string("legacy"))!="legacy")
    throw std::invalid_argument("self_energy_estimator currently supports imaginary-axis calculations only");
  if (e==Estimator::kugler && !kugler_available)
    throw std::invalid_argument("self_energy_estimator=kugler requires kugler=true and compatible analytic solver settings");
}
template<class T>
T choose(Estimator e,const T& dyson,const T& bulla,const T& kugler,double weight) {
  switch(e) {
  case Estimator::dyson: return dyson;
  case Estimator::bulla: return bulla;
  case Estimator::kugler: return kugler;
  case Estimator::blended: return (1.0-weight)*dyson+weight*bulla;
  }
  throw std::logic_error("invalid self-energy estimator");
}
}
#endif
