// Lightweight records shared by the legacy OCA enumerator and DLR evaluator.
//
// This deliberately contains no Eigen, Blitz, DLR, or integration state.  The
// expensive part of the legacy code is also the authoritative part: it walks
// the CIX operator graph, applies all selection rules, and evaluates the local
// matrix elements.  We keep that enumeration and export only its result in
// this small, stable format.  A compact evaluator can then reorganize the
// same paths without knowing anything about the CIX implementation.
#ifndef EDMFT_OCA_DLR_TYPES_H
#define EDMFT_OCA_DLR_TYPES_H

#include <array>

struct OCAPathRecord {
  // Pseudo-particle order around the Phi ring:
  //
  //      -----
  // p3 /  |    \ p0
  //   /   | s0  |
  //   |<--|---  |
  // p2\   |    / p1
  //    \  v s1/  
  //      -----
  //   p0 --s0--> p1 --s1--> p2 --s0--> p3 --s1--> p0 .
  //
  // The arrows only fix the operator/index convention; the two time-directed
  // hybridization lines are s0(t3-t1) and s1(t2).  This exact ordering is
  // used in every cut table in oca_dlr.h, so it should not be changed locally.
  std::array<int,4> pseudo{};

  // Independent hybridization-channel indices.  Signed/orbitally duplicated
  // output baths are reconstructed later, exactly as in legacy I_OCA.
  int s0=0, s1=0;

  // Local matrix-element products accumulated by the enumerator:
  //   0: mm               ordinary Phi / self-energy / physical-G weight
  //   1: nds1             F00 when hybridization s1 is cut
  //   2: nds0             F00 when hybridization s0 is cut
  //   3: nds1_anomalous   F01 when hybridization s1 is cut
  //   4: nds0_anomalous   F01 when hybridization s0 is cut
  // Several atomic paths can have identical pseudo and bath indices; Eigen's
  // sparse setFromTriplets later adds those duplicates exactly.
  std::array<double,5> weights{};
};

#endif
