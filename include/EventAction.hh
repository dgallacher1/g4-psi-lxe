#ifndef EventAction_hh
#define EventAction_hh 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <map>
#include <string>
#include <vector>

class G4Event;

// End of event: streams the LXe step hits, the fluence-plane hits, the per-event
// summary and (for geantino scans) the per-volume material budget.
//
// Ntuples (ids fixed in RunAction):
//   0 steps    1 events    2 fluence    3 budget    4 tree_chroma (one row per energy-depositing LXe step; GPU-simulation input)
//   5 tree_fuse (same rows + gamma interactions, lineage columns for Clustering.py)
//   6 primary_path (every step of the beam particle in every volume: the true path into the vessel)
class EventAction : public G4UserEventAction
{
 public:
  EventAction() = default;
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  // Called by SteppingAction (geantino scans only).
  void SetScanStart(G4double y, G4double z) { fY0 = y; fZ0 = z; }
  // Called by SteppingAction for trackID 1: one point of the primary's path (any volume).
  void AddPathPoint(const G4ThreeVector& p, G4double kinE, G4double t, const std::string& volume)
  { fPath.push_back({p.x(), p.y(), p.z(), kinE, t, volume}); }
  void AddBudget(const G4String& volume, G4double len, G4double x0, G4double lambda, G4double massThick);

 private:
  struct PathPt { G4double x, y, z, kinE, t; std::string vol; };
  std::vector<PathPt> fPath;
  struct Budget { G4double len = 0., x0 = 0., lambda = 0., massThick = 0.; };
  std::map<G4String, Budget> fBudget;
  G4double fY0 = 0., fZ0 = 0.;
  G4int fXenonHCID = -1;
  G4int fPlaneHCID = -1;
};

#endif
