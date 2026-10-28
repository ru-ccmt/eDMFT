// Per-stage defaults and exception-safe overrides for shared DLR workspaces.
#ifndef EDMFT_SUBBLOCK_POLICY_H
#define EDMFT_SUBBLOCK_POLICY_H
#include "json.hpp"
#include <stdexcept>
#include <string>
namespace edmft::imag_dlr {
struct SubblockPolicy {
  bool tca=true;
  bool fca=false;
};
inline SubblockPolicy resolve_subblocks(const nlohmann::json& value) {
  if(value.is_boolean()) {
    const bool enabled=value.get<bool>();
    return {enabled,enabled};
  }
  if(value.is_string() && value.get<std::string>()=="auto")return {};
  throw std::invalid_argument("higher_order_subblocks must be true, false, or auto");
}
// FCA MC normalization borrows TCA's accurate workspace. Restore its policy
// on all exit paths so a later TCA evaluation keeps its own setting.
template<class Workspace>
class ScopedAtomicOptimizations {
  Workspace& workspace;
  bool subblocks,skip_zero;
public:
  ScopedAtomicOptimizations(Workspace& w,bool enabled,bool zero)
    : workspace(w),subblocks(w.atomic_subblocks()),skip_zero(w.atomic_skip_zero_roots()) {
    workspace.set_atomic_optimizations(enabled,zero);
  }
  ScopedAtomicOptimizations(const ScopedAtomicOptimizations&)=delete;
  ScopedAtomicOptimizations& operator=(const ScopedAtomicOptimizations&)=delete;
  ~ScopedAtomicOptimizations() { workspace.set_atomic_optimizations(subblocks,skip_zero); }
};
}
#endif
