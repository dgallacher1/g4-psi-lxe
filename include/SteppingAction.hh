#ifndef SteppingAction_hh
#define SteppingAction_hh 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;
class G4GenericMessenger;
class G4ParticleDefinition;

// Material-budget scan.  Acts ONLY on geantinos (zero cost for physics runs
// beyond one pointer compare per step): accumulates, per logical volume,
// path length, X/X0, L/lambda_I and mass thickness, so the metal in the beam
// path can be checked against the CAD model.
//
//   /scan/stopAtLXe true|false   kill the geantino on entering the LXe (default true:
//                                 budget = everything the beam crosses BEFORE the xenon)
class SteppingAction : public G4UserSteppingAction
{
 public:
  explicit SteppingAction(EventAction* ea);
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step*) override;

 private:
  EventAction*                fEventAction;
  G4GenericMessenger*         fMessenger = nullptr;
  G4bool                      fStopAtLXe = true;
  const G4ParticleDefinition* fGeantino = nullptr;
  const G4ParticleDefinition* fChargedGeantino = nullptr;
};

#endif
