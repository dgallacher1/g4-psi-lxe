#ifndef PrimaryGeneratorAction_hh
#define PrimaryGeneratorAction_hh 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4GenericMessenger;
class G4Event;

// Two beam modes, switchable at runtime (no recompilation):
//
//   /beam/mode gps   G4GeneralParticleSource, configured in macros/beam.mac.
//   /beam/mode file  G4Beamline ASCII phase-space file (BLTrackFile), one primary
//                    per record.  The file is read by ONE shared, mutex-protected
//                    reader, so worker threads consume disjoint records.
//                    Use  /run/beamOn <number of records>.
//
// Frame mapping for BLTrackFile (G4Beamline: x horizontal, y vertical, z = beam):
//     sim.x = bl.z   (beam)       sim.y = bl.x      sim.z = bl.y   (vertical = vessel axis)
//     the same permutation is applied to the momentum, then
//     sim position += /beam/offsetX|Y|Z.
// This is a proper (right-handed) rotation.  Set the offset so the file's
// reference plane lies UPSTREAM of all material (see "beam/scan start plane
// suggestion" printed at startup).  If your vertical is not +z, say so and the
// mapping changes in GenerateFromFile().
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event*) override;

  // Master calls this at BeginOfRun so each run starts from record 0.
  static void RewindBeamFile();

 private:
  void DefineCommands();
  void GenerateFromFile(G4Event*);

  G4GeneralParticleSource* fGPS       = nullptr;
  G4GenericMessenger*      fMessenger = nullptr;

  G4String fMode     = "gps";
  G4String fBeamFile = "";
  G4double fOffX = 0., fOffY = 0., fOffZ = 0.;   // mm
};

#endif
