#include "RunAction.hh"

#include "G4AnalysisManager.hh"   // unified analysis (Geant4 >= 11)
#include "G4Run.hh"

RunAction::RunAction()
{
  auto* man = G4AnalysisManager::Instance();
  man->SetDefaultFileType("root");   // -> xenon_sim.root; 
  man->SetVerboseLevel(1);
  man->SetNtupleMerging(true);       // merge per-thread ntuples in MT mode

  // --- ntuple 0: steps (column ids run 0..13 in creation order) ---
  man->CreateNtuple("steps", "per-step path points and deposits in LXe");
  man->CreateNtupleIColumn("eventID");     // 0
  man->CreateNtupleIColumn("trackID");     // 1
  man->CreateNtupleIColumn("parentID");    // 2
  man->CreateNtupleIColumn("pdg");         // 3
  man->CreateNtupleDColumn("preX_mm");     // 4
  man->CreateNtupleDColumn("preY_mm");     // 5
  man->CreateNtupleDColumn("preZ_mm");     // 6
  man->CreateNtupleDColumn("postX_mm");    // 7
  man->CreateNtupleDColumn("postY_mm");    // 8
  man->CreateNtupleDColumn("postZ_mm");    // 9
  man->CreateNtupleDColumn("edep_keV");    // 10
  man->CreateNtupleDColumn("kinE_keV");    // 11
  man->CreateNtupleDColumn("stepLen_mm");  // 12
  man->CreateNtupleDColumn("time_ns");     // 13
  man->FinishNtuple();

  // --- ntuple 1: events ---
  man->CreateNtuple("events", "per-event totals in LXe");
  man->CreateNtupleIColumn("eventID");     // 0
  man->CreateNtupleDColumn("edep_keV");    // 1
  man->CreateNtupleIColumn("nSteps");      // 2
  man->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4AnalysisManager::Instance()->OpenFile("xenon_sim");
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}
