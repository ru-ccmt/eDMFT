// Canonical rooted topology definitions shared by every TCA implementation.
//
// Vertex numbers increase along the ordered pseudo-particle backbone.  Line
// zero is always the distinguished, unseparated hybridization attached to
// vertex zero.  Giving both topologies the same root (0,3) lets the DLR
// evaluator use one backbone algorithm; only the two separated chords differ.
#ifndef EDMFT_TCA_TOPOLOGY_H
#define EDMFT_TCA_TOPOLOGY_H

#include <array>
#include <stdexcept>

#include <blitz/array.h>

namespace edmft::tca {

using Endpoints=std::array<std::array<int,2>,3>;

inline constexpr std::array<Endpoints,2> topology_endpoints{{
  // A: cyclic relabeling of (0,2),(1,4),(3,5).
  {{{0,3},{2,4},{1,5}}},
  // B: the three mutually crossing diameter lines.
  {{{0,3},{1,4},{2,5}}}
}};

inline constexpr std::array<int,2> topology_crossing_sign{{+1,-1}};
// Rooting line (0,3) and fixing its direction removes the twofold arrow
// reversal counted by an undirected representation.  The fully crossing
// topology B nevertheless has three equivalent choices of the rooted
// diameter.  Its remaining automorphism therefore contributes 1/3.  This is
// also the factor obtained by comparison with the generic-diagram Monte Carlo
// evaluator; topology A has no residual multiplicity in this convention.
inline constexpr std::array<double,2> topology_symmetry_factor{{1.0,1.0/3.0}};
inline constexpr int root_line=0;

inline blitz::Array<int,2> MakeTopology(int topology)
{
  if (topology<0 || topology>=static_cast<int>(topology_endpoints.size()))
    throw std::out_of_range("TCA topology must be zero or one");
  blitz::Array<int,2> diag(3,2);
  for (int line=0;line<3;++line) {
    diag(line,0)=topology_endpoints[topology][line][0];
    diag(line,1)=topology_endpoints[topology][line][1];
  }
  return diag;
}

} // namespace edmft::tca

#endif
