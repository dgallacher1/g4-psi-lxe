#ifndef XenonHit_hh
#define XenonHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

// One hit == one step inside the LXe. Storing pre/post points preserves
// the full track polyline (the "path"); edep gives the deposit.
class XenonHit : public G4VHit
{
 public:
  XenonHit() = default;
  ~XenonHit() override = default;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  G4int    trackID  = -1;
  G4int    parentID = -1;
  G4int    pdg      = 0;
  G4double edep     = 0.;   // internal units (MeV)
  G4double kinE     = 0.;   // pre-step kinetic energy (MeV)
  G4double stepLen  = 0.;   // mm
  G4double gtime    = 0.;   // ns
  G4ThreeVector prePos;     // mm
  G4ThreeVector postPos;    // mm
};

using XenonHitsCollection = G4THitsCollection<XenonHit>;
extern G4ThreadLocal G4Allocator<XenonHit>* XenonHitAllocator;

inline void* XenonHit::operator new(size_t)
{
  if (!XenonHitAllocator) XenonHitAllocator = new G4Allocator<XenonHit>;
  return (void*)XenonHitAllocator->MallocSingle();
}

inline void XenonHit::operator delete(void* h)
{
  XenonHitAllocator->FreeSingle((XenonHit*)h);
}

#endif
