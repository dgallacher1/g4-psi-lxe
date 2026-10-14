#ifndef RunAction_hh
#define RunAction_hh 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class G4GenericMessenger;

// Creates the ntuples (defined in OutputSchema.hh) and owns the output file lifecycle and the
// /output/ commands (file tag and the switches of the optional ntuples).
class RunAction : public G4UserRunAction
{
 public:
  RunAction();
  ~RunAction() override;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  // /output/tag NAME : the next run writes xenon_sim_NAME.root (no dots in NAME).
  // Static so the worker threads see what the master's command set; it is only
  // written between runs.
  static G4String      sTag;
  G4GenericMessenger*  fMessenger = nullptr;
};

#endif
