#include "PrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include <sstream>

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fGPS = new G4GeneralParticleSource();
  DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if (fStream.is_open()) fStream.close();
  delete fGPS;
  delete fMessenger;
}

void PrimaryGeneratorAction::DefineCommands()
{
  fMessenger = new G4GenericMessenger(this, "/beam/", "Primary beam control");
  fMessenger->DeclareProperty("mode", fMode,
      "Beam mode: 'gps' (default) or 'file' (G4Beamline BLTrackFile)");
  fMessenger->DeclareProperty("file", fBeamFile,
      "Path to a G4Beamline ASCII phase-space file (used when mode=file)");
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt)
{
  if (fMode == "file")
    GenerateFromFile(evt);
  else
    fGPS->GeneratePrimaryVertex(evt);
}

// G4Beamline BLTrackFile columns:
//   x y z Px Py Pz t PDGid EventID TrackID ParentID Weight
// Units: mm, MeV/c, ns. Header/comment lines start with '#'.
void PrimaryGeneratorAction::GenerateFromFile(G4Event* evt)
{
  if (!fStreamOpen) {
    fStream.open(fBeamFile);
    if (!fStream.is_open()) {
      G4Exception("PrimaryGeneratorAction::GenerateFromFile", "beam_file_open",
                  FatalException, ("cannot open " + fBeamFile).c_str());
      return;
    }
    fStreamOpen = true;
  }

  std::string line;
  while (std::getline(fStream, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    G4double x, y, z, px, py, pz, t, w;
    G4int pdg, evid, tid, pid;
    if (!(iss >> x >> y >> z >> px >> py >> pz >> t >> pdg >> evid >> tid >> pid >> w))
      continue;

    auto* def = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
    if (!def) continue;

    auto* vtx = new G4PrimaryVertex(G4ThreeVector(x * mm, y * mm, z * mm), t * ns);
    vtx->SetPrimary(new G4PrimaryParticle(def, px * MeV, py * MeV, pz * MeV));
    evt->AddPrimaryVertex(vtx);
    return;  // one line -> one primary -> one event
  }

  // Ran out of phase-space records: stop cleanly.
  G4RunManager::GetRunManager()->AbortRun();
}
