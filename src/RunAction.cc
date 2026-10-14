#include "RunAction.hh"
#include "GeometryParameters.hh"
#include "PrimaryGeneratorAction.hh"
#include "OutputOptions.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4Threading.hh"
#include "G4GenericMessenger.hh"

G4String RunAction::sTag = "";

RunAction::~RunAction() { delete fMessenger; }

RunAction::RunAction()
{
  if (G4Threading::IsMasterThread()) {      // master only: the commands are not broadcast
    fMessenger = new G4GenericMessenger(this, "/output/", "Output file control");
    fMessenger->DeclareProperty("tag", sTag,
        "Suffix for the next run's ROOT file: xenon_sim_<tag>.root (no dots)")
        .SetToBeBroadcasted(false);
    fMessenger->DeclareProperty("chroma", gWriteChroma,
        "false: do not fill tree_chroma (GPU-simulation input)").SetToBeBroadcasted(false);
    fMessenger->DeclareProperty("fuse", gWriteFuse,
        "false: do not fill tree_fuse (lineage table read by Clustering.load_interactions)")
        .SetToBeBroadcasted(false);
    fMessenger->DeclareProperty("steps", gWriteSteps,
        "false: do not fill the steps ntuple (smaller files when only tree_chroma is needed)")
        .SetToBeBroadcasted(false);
  }

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
  man->CreateNtupleDColumn("niel_keV");    // 14  non-ionising part of edep_keV (nuclear stopping)
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
  man->CreateNtupleDColumn("vtxX_mm");     // 13  track vertex (where the particle was created)
  man->CreateNtupleDColumn("vtxY_mm");     // 14
  man->CreateNtupleDColumn("vtxZ_mm");     // 15
  man->CreateNtupleSColumn("creator");     // 16  creator process, "primary" for the beam particle
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

  // --- ntuple 4: tree_chroma, the input of the GPU optical-photon simulation ----------------
  // Same layout as the tree_chroma of the other Geant4 simulation that feeds it: one ROW PER
  // ENERGY-DEPOSITING LXe STEP, every column a double. Rows of one event share "event"; "hits" is
  // the index of the row within the event. Units: mm, ns, keV (edep, energy, momentum).
  //   particle_name is the numeric class used there: 1 gamma, 2 neutron, 3 e-, 4 e+, 5 other
  //   (the real identity is in PDG).
  man->CreateNtuple("tree_chroma", "Scintillation each hit info");
  for (const char* c : {"event", "hits", "xpos", "ypos", "zpos", "time", "edep", "energy",
                        "x_momentum", "y_momentum", "z_momentum", "particle_name", "PDG", "trackID"})
    man->CreateNtupleDColumn(c);                                 // columns 0..13 in this order
  man->FinishNtuple();

  // --- ntuple 5: tree_fuse, the columns Clustering.load_interactions reads ------------------
  // One row per energy-depositing LXe step, plus the zero-deposit gamma interaction steps (they are
  // the lineage parents of the electrons they release). Units: keV, ns, mm.
  man->CreateNtuple("tree_fuse", "per-step lineage table for Clustering.py");
  man->CreateNtupleIColumn("eventid");       // 0
  man->CreateNtupleIColumn("trackid");       // 1
  man->CreateNtupleIColumn("parentid");      // 2
  man->CreateNtupleSColumn("type");          // 3  particle name ("e-", "gamma", "alpha", "Pb208[...]")
  man->CreateNtupleSColumn("parenttype");    // 4  parent's particle name, "none" for primaries
  man->CreateNtupleSColumn("creaproc");      // 5  creator process, "none" for primaries
  man->CreateNtupleSColumn("edproc");        // 6  process that ended the step
  man->CreateNtupleDColumn("ed");            // 7  keV
  man->CreateNtupleDColumn("t");             // 8  ns (post-step)
  man->CreateNtupleDColumn("xp");            // 9  mm (post-step)
  man->CreateNtupleDColumn("yp");            // 10
  man->CreateNtupleDColumn("zp");            // 11
  man->CreateNtupleDColumn("xp_pri");        // 12 primary vertex, mm
  man->CreateNtupleDColumn("yp_pri");        // 13
  man->CreateNtupleDColumn("zp_pri");        // 14
  // Same values again under the plain names: Clustering.load_interactions asks uproot for all of
  // xp/x, xp_pri/x_pri... in one call and uproot raises on a missing branch; it also only keeps
  // x_pri when the plain-named column exists. Prefers xp over x, so no double counting.
  for (const char* c : {"x", "y", "z", "x_pri", "y_pri", "z_pri"})
    man->CreateNtupleDColumn(c);             // 15..20, mm
  man->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
  // Master only: restart a shared G4Beamline phase-space file so every run
  // consumes it from the first record.
  if (G4Threading::IsMasterThread()) PrimaryGeneratorAction::RewindBeamFile();

  // Placeholder geometry can never be mistaken for CAD output: the file name says so.
  G4String name = Geo::isDemo ? "xenon_sim_DEMO" : "xenon_sim";
  if (!sTag.empty()) name += "_" + sTag;
  G4AnalysisManager::Instance()->OpenFile(name);
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto* man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();
}
