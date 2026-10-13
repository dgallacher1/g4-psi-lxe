#ifndef XenonSD_hh
#define XenonSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "XenonHit.hh"

class G4Step;
class G4HCofThisEvent;

class XenonSD : public G4VSensitiveDetector
{
 public:
  explicit XenonSD(const G4String& name);
  ~XenonSD() override = default;

  void   Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

 private:
  XenonHitsCollection* fHC   = nullptr;
  G4int                fHCID = -1;
};

#endif
