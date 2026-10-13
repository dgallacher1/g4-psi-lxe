#include "PrimaryGeneratorAction.hh"

#include "G4GeneralParticleSource.hh"
#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4AutoLock.hh"
#include "G4SystemOfUnits.hh"

#include <atomic>
#include <fstream>
#include <sstream>

namespace {

// One reader shared by all threads.
G4Mutex       gFileMutex = G4MUTEX_INITIALIZER;
std::ifstream gStream;
G4String      gOpenName;
long          gRecordsRead = 0;
std::atomic<bool> gWarnedEOF{false};

}  // namespace

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fGPS = new G4GeneralParticleSource();
  DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fGPS;
  delete fMessenger;
}

void PrimaryGeneratorAction::RewindBeamFile()
{
  G4AutoLock lock(&gFileMutex);
  if (gStream.is_open()) gStream.close();
  gStream.clear();
  gOpenName.clear();
  gRecordsRead = 0;
  gWarnedEOF = false;
}

void PrimaryGeneratorAction::DefineCommands()
{
  fMessenger = new G4GenericMessenger(this, "/beam/", "Primary beam control");
  fMessenger->DeclareProperty("mode", fMode,
      "Beam mode: 'gps' (default) or 'file' (G4Beamline BLTrackFile)");
  fMessenger->DeclareProperty("file", fBeamFile,
      "Path to a G4Beamline ASCII phase-space file (used when mode=file)");
  fMessenger->DeclarePropertyWithUnit("offsetX", "mm", fOffX,
      "Sim-frame x of the BLTrackFile origin (file z=0 plane)");
  fMessenger->DeclarePropertyWithUnit("offsetY", "mm", fOffY, "Sim-frame y of the BLTrackFile origin");
  fMessenger->DeclarePropertyWithUnit("offsetZ", "mm", fOffZ, "Sim-frame z of the BLTrackFile origin");
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* evt)
{
  if (fMode == "file")
    GenerateFromFile(evt);
  else
    fGPS->GeneratePrimaryVertex(evt);
}

// G4Beamline BLTrackFile columns:
//   x y z Px Py Pz t PDGid EventID TrackID ParentID Weight      [mm, MeV/c, ns]
void PrimaryGeneratorAction::GenerateFromFile(G4Event* evt)
{
  std::string line;
  G4double x = 0, y = 0, z = 0, px = 0, py = 0, pz = 0, t = 0, w = 0;
  G4int pdg = 0, evid = 0, tid = 0, pid = 0;
  G4ParticleDefinition* def = nullptr;
  G4bool got = false;

  {
    G4AutoLock lock(&gFileMutex);
    if (!gStream.is_open() || gOpenName != fBeamFile) {
      if (gStream.is_open()) gStream.close();
      gStream.clear();
      gStream.open(fBeamFile);
      if (!gStream.is_open())
        G4Exception("PrimaryGeneratorAction::GenerateFromFile", "beam_file_open", FatalException,
                    ("cannot open beam file '" + fBeamFile + "' (set /beam/file)").c_str());
      gOpenName = fBeamFile;
      gRecordsRead = 0;
    }
    while (std::getline(gStream, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream iss(line);
      if (!(iss >> x >> y >> z >> px >> py >> pz >> t >> pdg >> evid >> tid >> pid >> w)) continue;
      def = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
      if (!def) continue;
      ++gRecordsRead;
      got = true;
      break;
    }
  }

  if (!got) {
    // Out of records: emit no primary (empty event). Warn once.
    if (!gWarnedEOF.exchange(true))
      G4Exception("PrimaryGeneratorAction::GenerateFromFile", "beam_file_eof", JustWarning,
                  "Beam file exhausted: remaining events are empty. Use /run/beamOn <number of records>.");
    return;
  }

  // BLTrackFile (x,y,z) -> sim (x,y,z) = (bl.z, bl.x, bl.y), + offset.
  const G4ThreeVector pos(z * mm + fOffX, x * mm + fOffY, y * mm + fOffZ);
  const G4ThreeVector mom(pz * MeV, px * MeV, py * MeV);

  auto* vtx = new G4PrimaryVertex(pos, t * ns);
  vtx->SetPrimary(new G4PrimaryParticle(def, mom.x(), mom.y(), mom.z()));
  evt->AddPrimaryVertex(vtx);
}
