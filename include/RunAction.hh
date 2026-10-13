#ifndef RunAction_hh
#define RunAction_hh 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

// Owns the g4tools/G4AnalysisManager ntuple schema and file lifecycle.
// Two ntuples:
//   0 "steps"  : one row per step in the LXe (path points + edep)
//   1 "events" : one row per event (total edep, step multiplicity)
class RunAction : public G4UserRunAction
{
 public:
  RunAction();
  ~RunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;
};

#endif
