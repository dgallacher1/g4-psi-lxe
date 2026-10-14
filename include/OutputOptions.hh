#ifndef OutputOptions_hh
#define OutputOptions_hh 1

#include "globals.hh"

// Output switches, set on the master between runs (/output/chroma, /output/steps) and only read
// during a run, so one shared copy is enough.
inline G4bool gWriteChroma = true;    // fill the tree_chroma ntuple (input of the GPU simulation)
inline G4bool gWriteSteps  = true;    // fill the per-step "steps" ntuple

#endif
