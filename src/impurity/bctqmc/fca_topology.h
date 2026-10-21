// Fourth-order skeletons, in the user's cyclic vertex convention.
// IMPORTANT: these weights are for ALL 16 arrow assignments. No arrow is
// fixed and no pseudo/physical cut is discarded in the reference evaluator.
#ifndef EDMFT_FCA_TOPOLOGY_H
#define EDMFT_FCA_TOPOLOGY_H
#include <array>
#include <stdexcept>
#include <algorithm>
#include <vector>

namespace edmft::fca {
using Endpoints=std::array<std::array<int,2>,4>;
inline constexpr std::array<Endpoints,6> endpoints{{
  {{{0,4},{1,7},{2,6},{3,5}}}, // A
  {{{0,4},{1,7},{2,5},{3,6}}}, // B
  {{{0,4},{1,6},{2,5},{3,7}}}, // C
  {{{0,4},{1,5},{2,6},{3,7}}}, // D: repeated (2,6) removed
  {{{0,3},{1,6},{2,4},{5,7}}}, // E
  {{{0,3},{1,6},{2,5},{4,7}}}  // F
}};
inline constexpr std::array<double,6> symmetry{{0.5,1.0,0.5,0.125,1.0,0.25}};
// Chord crossing parity only. Atomic-operator signs and signed bath/Sindex
// labels are handled separately, exactly as in the generic MC/TCA code.
inline constexpr std::array<int,6> crossing_sign{{-1,+1,-1,+1,-1,+1}};
inline constexpr const char* names="ABCDEF";

struct Line { int original=-1,a=-1,b=-1,d=0; };
struct Cut {
  Line root;
  std::array<Line,3> split;
  int q=0;
};

// Rotate a pseudo cut (edge c -> c+1) to the closing edge. The integration
// root is the line incident at the NEW origin, not necessarily original line
// zero. Endpoint reversal flips an arrow bit, not its bath label/residue.
inline Cut self_cut(int topology,int cut,const std::array<int,4>& direction)
{
  if (topology<0 || topology>=6 || cut<0 || cut>=8)
    throw std::out_of_range("invalid FCA pseudo cut");
  Cut result;
  const int start=(cut+1)%8;
  int n=0;
  for (int l=0;l<4;++l) {
    int a=(endpoints[topology][l][0]-start+8)%8;
    int b=(endpoints[topology][l][1]-start+8)%8;
    int d=direction[l];
    if (d<0 || d>1) throw std::out_of_range("invalid FCA direction");
    if (a>b) { std::swap(a,b); d^=1; }
    Line line{l,a,b,d};
    if (a==0) { result.root=line; result.q=b; }
    else result.split[n++]=line;
  }
  if (n!=3 || result.q==0) throw std::logic_error("invalid FCA matching");
  return result;
}

// Rotate a physical cut so the creation boundary is zero, matching TCA's
// ordinary/Bulla contraction. There is NO Delta on this external line.
inline Cut physical_cut(int topology,int cut,const std::array<int,4>& direction)
{
  if (topology<0 || topology>=6 || cut<0 || cut>=4)
    throw std::out_of_range("invalid FCA physical cut");
  const auto edge=endpoints[topology][cut];
  const int start=direction[cut]==0 ? edge[1] : edge[0];
  const int finish=direction[cut]==0 ? edge[0] : edge[1];
  Cut result;
  result.q=(finish-start+8)%8;
  int n=0;
  for (int l=0;l<4;++l) if (l!=cut) {
    int a=(endpoints[topology][l][0]-start+8)%8;
    int b=(endpoints[topology][l][1]-start+8)%8;
    int d=direction[l];
    if (d<0 || d>1) throw std::out_of_range("invalid FCA direction");
    if (a>b) { std::swap(a,b); d^=1; }
    result.split[n++]={l,a,b,d};
  }
  return result;
}

// Two rooted cuts are interchangeable only when their ordered-time chord
// geometry is identical after rotation and relabeling the three summed
// internal lines. This is NOT a comparison of symmetrized pseudo labels.
// All bath, arrow, and pole assignments must still be summed for each
// representative. In particular, one fixed pole triple need not be equal
// before the internal-line labels have also been permuted.
struct CutShape {
  int q;
  std::array<std::array<int,2>,3> internal;
  bool operator==(const CutShape&) const = default;
};
inline CutShape shape(const Cut& cut) {
  CutShape key; key.q=cut.q;
  for (int j=0;j<3;++j) key.internal[j]={{cut.split[j].a,cut.split[j].b}};
  std::sort(key.internal.begin(),key.internal.end());
  return key;
}
struct CutOrbit {
  int topology=0, representative=0;
  std::vector<int> members;
};
inline std::vector<CutOrbit> cut_orbits(bool physical,bool reduced) {
  std::vector<CutOrbit> result;
  for (int t=0;t<6;++t) {
    std::vector<CutShape> keys;
    std::vector<std::size_t> positions;
    for (int c=0;c<8;++c) {
      // Physical index encodes (cut line, orientation of external line).
      // Every external endpoint is kept; the creation boundary is always
      // rotated to zero by physical_cut(), preserving the Bulla insertion.
      std::array<int,4> arrows{};
      if (physical) arrows[c/2]=c%2;
      const auto key=shape(physical ? physical_cut(t,c/2,arrows)
                                    : self_cut(t,c,arrows));
      const auto found=reduced ? std::find(keys.begin(),keys.end(),key) : keys.end();
      if (found!=keys.end()) {
        result[positions[found-keys.begin()]].members.push_back(c);
      } else {
        keys.push_back(key); positions.push_back(result.size());
        result.push_back({t,c,{c}});
      }
    }
  }
  return result;
}
} // namespace edmft::fca
#endif
