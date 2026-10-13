#ifndef PlaneHit_hh
#define PlaneHit_hh 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "globals.hh"

// One charged-particle step inside a fluence plane.
//
// Fluence can be formed offline two ways from these rows:
//   * crossing estimator : sum over rows with entering==1 of 1/|cosTheta|, / area
//   * track-length       : sum(stepLen_mm) / (area * thickness)
// (they agree for a thin plane; the second is exact and also counts particles
// born inside the plane).
class PlaneHit : public G4VHit
{
 public:
  PlaneHit() = default;
  ~PlaneHit() override = default;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  G4int    planeID  = -1;
  G4int    trackID  = -1;
  G4int    parentID = -1;
  G4int    pdg      = 0;
  G4int    entering = 0;     // 1 if the step starts on the plane's boundary (i.e. a crossing)
  G4double u        = 0.;    // local in-plane coordinates of the pre-step point [mm]
  G4double v        = 0.;
  G4double kinE     = 0.;    // pre-step kinetic energy
  G4double cosTheta = 0.;    // pre-step direction . plane normal (signed)
  G4double stepLen  = 0.;
  G4double edep     = 0.;
  G4double gtime    = 0.;
};

using PlaneHitsCollection = G4THitsCollection<PlaneHit>;
extern G4ThreadLocal G4Allocator<PlaneHit>* PlaneHitAllocator;

inline void* PlaneHit::operator new(size_t)
{
  if (!PlaneHitAllocator) PlaneHitAllocator = new G4Allocator<PlaneHit>;
  return (void*)PlaneHitAllocator->MallocSingle();
}

inline void PlaneHit::operator delete(void* h)
{
  PlaneHitAllocator->FreeSingle((PlaneHit*)h);
}

#endif
