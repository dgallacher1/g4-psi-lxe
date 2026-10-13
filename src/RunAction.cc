#include "RunAction.hh"
#include "GeometryParameters.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4Threading.hh"

RunAction::RunAction()
{
  auto* man = G4AnalysisManager::Instance();
  man->SetDefaultFileType("root");
  man->SetVerboseLevel(1);
  man->SetNtupleMerging(true);       // merge per-thread ntuples in MT mode

  // --- ntuple 0: steps in the LXe (column ids = creation order) ---
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
  man->CreateNtuple("events", "per-event totals");
  man->CreateNtupleIColumn("eventID");          // 0
  man->CreateNtupleDColumn("edep_keV");         // 1  bulk LXe only
  man->CreateNtupleIColumn("nSteps");           // 2
  man->CreateNtupleDColumn("edepPlanes_keV");   // 3  deposits inside the fluence planes
  man->FinishNtuple();

  // --- ntuple 2: charged-particle steps in the fluence planes ---
  man->CreateNtuple("fluence", "charged-particle steps in the virtual SiPM planes");
  man->CreateNtupleIColumn("eventID");     // 0
  man->CreateNtupleIColumn("planeID");     // 1  index into Geo::planes / geometry_manifest.json
  man->CreateNtupleIColumn("trackID");     // 2
  man->CreateNtupleIColumn("parentID");    // 3
  man->CreateNtupleIColumn("pdg");         // 4
  man->CreateNtupleIColumn("entering");    // 5  1 = step starts on the plane boundary
  man->CreateNtupleDColumn("u_mm");        // 6  local in-plane coordinates
  man->CreateNtupleDColumn("v_mm");        // 7
  man->CreateNtupleDColumn("kinE_keV");    // 8  pre-step
  man->CreateNtupleDColumn("cosTheta");    // 9  direction . plane normal
  man->CreateNtupleDColumn("stepLen_mm");  // 10
  man->CreateNtupleDColumn("edep_keV");    // 11
  man->CreateNtupleDColumn("time_ns");     // 12
  man->FinishNtuple();

  // --- ntuple 3: geantino material budget, one row per (event, volume) ---
  man->CreateNtuple("budget", "geantino material budget per logical volume");
  man->CreateNtupleIColumn("eventID");        // 0
  man->CreateNtupleSColumn("volume");         // 1
  man->CreateNtupleDColumn("y0_mm");          // 2  geantino start point
  man->CreateNtupleDColumn("z0_mm");          // 3
  man->CreateNtupleDColumn("length_mm");      // 4
  man->CreateNtupleDColumn("X0");             // 5  path length / radiation length
  man->CreateNtupleDColumn("lambdaI");        // 6  path length / nuclear interaction length
  man->CreateNtupleDColumn("massThick_gcm2"); // 7
  man->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  // Master only: restart a shared G4Beamline phase-space file so every run
  // consumes it from the first record.
  if (G4Threading::IsMasterThread()) PrimaryGeneratorAction::RewindBeamFile();

  // Placeholder geometry can never be mistaken for CAD output: the file name says so.
  G4AnalysisManager::Instance()->OpenFile(Geo::isDemo ? "xenon_sim_DEMO" : "xenon_sim");
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}
