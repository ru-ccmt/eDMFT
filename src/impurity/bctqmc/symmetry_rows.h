#ifndef EDMFT_SYMMETRY_ROWS_H
#define EDMFT_SYMMETRY_ROWS_H
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <utility>
namespace edmft::time_reversal {
// Rows are shortened only AFTER the unique two-sided root insertion. Every
// subsequent operation multiplies on the right. This is an output projector,
// not an approximation to a propagator or an internal atomic state space.
struct Rows {
  std::vector<std::vector<int>> index;
  std::vector<std::vector<int>> weight;
  int full=0,retained=0;
  bool active() const {return retained>0 && retained<full;}
  template<class Function> void project(Function& function) const {
    if(!active())return;
    for(auto it=function.begin();it!=function.end();) {
      const auto& rows=index.at(it->first.first);
      if(rows.empty()){it=function.erase(it);continue;}
      for(auto& values:it->second) {
        std::decay_t<decltype(values)> selected(rows.size(),values.cols());
        for(std::size_t i=0;i<rows.size();++i)selected.row(i)=values.row(rows[i]);
        values=std::move(selected);
      }
      ++it;
    }
  }
};
}
#endif
