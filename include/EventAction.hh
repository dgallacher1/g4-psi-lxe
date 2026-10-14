#ifndef EventAction_hh
#define EventAction_hh 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <map>
#include <string>
#include <vector>

class G4Event;

// End of event: writes the event's outputs. The ntuples and columns are defined in OutputSchema.hh;
// each output has its own Fill... function below, gated by its /output/ switch.
class G4HCofThisEvent;
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
  struct Totals { G4double lxeEdep = 0., planeEdep = 0.; G4int nSteps = 0; };

  // hit collections of the event -> steps, tree_chroma, tree_fuse (one pass over the LXe hits)
  void FillLXeHits(G4HCofThisEvent*, G4int eid, const G4ThreeVector& primaryVertex, Totals&);
  void FillPlaneHits(G4HCofThisEvent*, G4int eid, Totals&);
  void FillPrimaryPath(G4int eid) const;                 // primary_path
  void FillEventSummary(G4int eid, const Totals&) const; // events
  void FillBudget(G4int eid) const;                      // budget (geantino scans)

  struct PathPt { G4double x, y, z, kinE, t; std::string vol; };
  std::vector<PathPt> fPath;
  struct Budget { G4double len = 0., x0 = 0., lambda = 0., massThick = 0.; };
  std::map<G4String, Budget> fBudget;
  G4double fY0 = 0., fZ0 = 0.;
  G4int fXenonHCID = -1;
  G4int fPlaneHCID = -1;
};

#endif
