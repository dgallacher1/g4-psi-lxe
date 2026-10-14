#ifndef XenonHit_hh
#define XenonHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "TrackLabels.hh"
#include <string>

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
  G4double edep     = 0.;   // total deposit, internal units (MeV); includes niel
  G4double niel     = 0.;   // non-ionising part of edep (MeV): nuclear stopping / recoils
  G4double kinE     = 0.;   // pre-step kinetic energy (MeV)
  G4double stepLen  = 0.;   // mm
  G4double gtime    = 0.;   // ns, pre-step
  G4double postTime = 0.;   // ns, post-step global time
  G4ThreeVector prePos;     // mm
  G4ThreeVector postPos;    // mm
  G4ThreeVector preMom;     // pre-step momentum (MeV/c)
  const TrackLabel* label = nullptr;   // type / parent type / creator process of the track (this event)
  std::string edproc;       // process that ended the step ("compt", "eIoni", "Transportation" ...)
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
