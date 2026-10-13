#ifndef PlaneSD_hh
#define PlaneSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "G4RotationMatrix.hh"
#include "PlaneHit.hh"
#include "GeometryParameters.hh"

#include <array>

class G4Step;
class G4HCofThisEvent;

// Records every step of a CHARGED particle inside any fluence plane.
// The plane is identified by the placement copy number; local coordinates and
// the plane normal are rebuilt from Geo::planes[] so the SD and the geometry
// can never disagree about orientation.
class PlaneSD : public G4VSensitiveDetector
{
 public:
  explicit PlaneSD(const G4String& name);
  ~PlaneSD() override = default;

  void   Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

 private:
  PlaneHitsCollection* fHC   = nullptr;
  G4int                fHCID = -1;
  std::array<G4RotationMatrix, Geo::nPlanes> fRot;
};

#endif
