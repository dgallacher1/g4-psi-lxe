#ifndef PhysicsList_hh
#define PhysicsList_hh 1

#include "G4VModularPhysicsList.hh"

// Modular EM-centric list. Deliberately NO optical physics
// (no G4OpticalPhysics) — we transport and score charged/neutral
// particles and energy deposits only, per the requirements.
class PhysicsList : public G4VModularPhysicsList
{
 public:
  PhysicsList();
  ~PhysicsList() override = default;
};

#endif
