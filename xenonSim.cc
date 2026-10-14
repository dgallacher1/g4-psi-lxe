// xenonSim.cc — application entry point.
//
// Batch:        ./xenonSim run.mac
// Interactive:  ./xenonSim            (opens the vis viewer via vis.mac)

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "Materials.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include <ctime>

int main(int argc, char** argv)
{
  // Interactive session only when no macro is supplied.
  G4UIExecutive* ui = (argc == 1) ? new G4UIExecutive(argc, argv) : nullptr;

  // Seed from wall clock; override in a macro with /random/setSeeds for reproducibility.
  CLHEP::HepRandom::setTheSeed(static_cast<long>(time(nullptr)));

  // Default = MT if Geant4 was built multithreaded, otherwise Serial.
  auto* runManager =
      G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

  Materials::Define();                        // SS304, TiWindowAlloy: must exist before the geometry is built
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization());

  auto* visManager = new G4VisExecutive;
  visManager->Initialize();

  auto* UI = G4UImanager::GetUIpointer();
  UI->SetMaxHistSize(5000);                   // keep the executed commands (default 20) for run_info

  if (ui) {                                   // interactive
    UI->ApplyCommand("/control/execute macros/vis.mac");
    ui->SessionStart();
    delete ui;
  } else {                                    // batch
    UI->ApplyCommand("/control/execute " + G4String(argv[1]));
  }

  delete visManager;
  delete runManager;
  return 0;
}
