#include "RunAction.hh"
#include "GeometryParameters.hh"
#include "PrimaryGeneratorAction.hh"
#include "OutputOptions.hh"
#include "OutputSchema.hh"
#include "RunMeta.hh"

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
    fMessenger->DeclareProperty("path", gWritePath,
        "false: do not fill primary_path (beam particle steps in every volume)").SetToBeBroadcasted(false);
    fMessenger->DeclareProperty("steps", gWriteSteps,
        "false: do not fill the steps ntuple (smaller files when only tree_chroma is needed)")
        .SetToBeBroadcasted(false);
  }

  auto* man = G4AnalysisManager::Instance();
  man->SetDefaultFileType("root");
  man->SetVerboseLevel(1);
  man->SetNtupleMerging(true);       // merge per-thread ntuples in MT mode

  // The ntuples and their columns are defined in OutputSchema.hh. Create them in id order.
  Out::Create<Out::Steps>(man);
  Out::Create<Out::Events>(man);
  Out::Create<Out::Fluence>(man);
  Out::Create<Out::Budget>(man);
  Out::Create<Out::Chroma>(man);    // input of the GPU optical-photon simulation
  Out::Create<Out::Fuse>(man);      // input of Clustering.py (lineage clustering)
  Out::Create<Out::Path>(man);      // beam particle path through every volume
  Out::Create<Out::Meta>(man);      // provenance of the run (RunMeta)
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  // Master only: restart a shared G4Beamline phase-space file so every run
  // consumes it from the first record.
  if (G4Threading::IsMasterThread()) PrimaryGeneratorAction::RewindBeamFile();

  // Placeholder geometry can never be mistaken for CAD output: the file name says so.
  G4String name = Geo::isDemo ? "xenon_sim_DEMO" : "xenon_sim";
  if (!sTag.empty()) name += "_" + sTag;
  // Provenance: collected once on the master, written to <name>_run_info.json at once.
  if (G4Threading::IsMasterThread()) RunMeta::Collect(run, name);
  G4AnalysisManager::Instance()->OpenFile(name);
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto* man = G4AnalysisManager::Instance();

  // run_info: one thread fills it (worker 0 in MT, the only thread otherwise); the merge on Write()
  // brings the rows into the master file. Rows come from RunMeta, collected by the master at the start.
  if (!G4Threading::IsMultithreadedApplication() || G4Threading::G4GetThreadId() == 0) {
    Out::Row<Out::Meta> r(man);
    for (const auto& kv : RunMeta::Get()) {
      r.set(Out::Meta::key, kv.first);
      r.set(Out::Meta::value, kv.second);
      r.add();
    }
  }
  man->Write();
  man->CloseFile();
}
