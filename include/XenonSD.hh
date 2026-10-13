#ifndef XenonSD_hh
#define XenonSD_hh 1

#include "G4VSensitiveDetector.hh"
#include "XenonHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4GenericMessenger;

class XenonSD : public G4VSensitiveDetector
{
 public:
  explicit XenonSD(const G4String& name);
  ~XenonSD() override;

  void   Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override;

 private:
  XenonHitsCollection* fHC   = nullptr;
  G4int                fHCID = -1;

  // /lxe/onlyEdep true : store only steps with edep > 0 (smaller files; fine for a
  // photon-yield consumer). Default false keeps every step, so the neutral fly-throughs
  // and the full trajectory polyline that lxe_look.py draws are there.
  G4bool               fOnlyEdep = false;
  G4GenericMessenger*  fMessenger = nullptr;
};

#endif
