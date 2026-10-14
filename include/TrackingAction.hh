#ifndef TrackingAction_hh
#define TrackingAction_hh 1

#include "G4UserTrackingAction.hh"

// Records, for every track, its particle name, its parent's particle name and its creator
// process (TrackLabels.hh). The "tree_fuse" ntuple needs them for the lineage clustering.
class TrackingAction : public G4UserTrackingAction
{
 public:
  void PreUserTrackingAction(const G4Track*) override;
};

#endif
