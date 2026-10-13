#ifndef ActionInitialization_hh
#define ActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"

class ActionInitialization : public G4VUserActionInitialization
{
 public:
  ActionInitialization() = default;
  ~ActionInitialization() override = default;

  void BuildForMaster() const override;  // master thread: RunAction only
  void Build() const override;           // worker threads: full set
};

#endif
