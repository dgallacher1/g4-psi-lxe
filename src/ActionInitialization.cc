#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction());
}

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGeneratorAction());
  SetUserAction(new RunAction());
  auto* ea = new EventAction();
  SetUserAction(ea);
  SetUserAction(new SteppingAction(ea));   // acts on geantinos only
  SetUserAction(new TrackingAction());     // particle / parent / creator labels for tree_fuse
}
