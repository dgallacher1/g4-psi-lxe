#ifndef TrackLabels_hh
#define TrackLabels_hh 1

#include "globals.hh"
#include <string>
#include <unordered_map>

// Particle labels per track of the current event, filled by TrackingAction when a track starts
// and read by XenonSD / EventAction (the lineage columns of the "tree_fuse" ntuple).
// Thread-local: every worker thread tracks its own events. Cleared at the start of each event;
// the map is node-based, so pointers into it stay valid for the whole event.
struct TrackLabel
{
  std::string type;        // particle name ("e-", "gamma", "alpha", "Np237[59.540]" ...)
  std::string parenttype;  // name of the parent track's particle, "none" for a primary
  std::string creaproc;    // creator process, "none" for a primary
};

inline thread_local std::unordered_map<G4int, TrackLabel> gTrackLabels;

#endif
