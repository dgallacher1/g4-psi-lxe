#include "SteppingAction.hh"
#include "EventAction.hh"
#include "OutputOptions.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Geantino.hh"
#include "G4ChargedGeantino.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(EventAction* ea) : fEventAction(ea)
{
  fGeantino        = G4Geantino::Definition();
  fChargedGeantino = G4ChargedGeantino::Definition();

  fMessenger = new G4GenericMessenger(this, "/scan/", "Material-budget scan control");
  fMessenger->DeclareProperty("stopAtLXe", fStopAtLXe,
      "Kill geantinos on entering the LXe (budget up to the xenon)");
}

SteppingAction::~SteppingAction() { delete fMessenger; }

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  auto* track = step->GetTrack();
  const auto* def = track->GetDefinition();

  // Beam particle path in every volume (air, foils, vacuum, LXe ...): start point, then the post-step
  // point of every step labelled with the volume the step ran in. Written to the "primary_path" ntuple.
  if (gWritePath && track->GetTrackID() == 1 && def != fGeantino && def != fChargedGeantino) {
    auto* pv0 = step->GetPreStepPoint()->GetPhysicalVolume();
    const std::string vol = pv0 ? pv0->GetName() : "OutOfWorld";
    if (track->GetCurrentStepNumber() == 1)
      fEventAction->AddPathPoint(step->GetPreStepPoint()->GetPosition(), step->GetPreStepPoint()->GetKineticEnergy(),
                                 step->GetPreStepPoint()->GetGlobalTime(), vol);
    fEventAction->AddPathPoint(step->GetPostStepPoint()->GetPosition(), step->GetPostStepPoint()->GetKineticEnergy(),
                               step->GetPostStepPoint()->GetGlobalTime(), vol);
  }

  if (def != fGeantino && def != fChargedGeantino) return;

  auto* pre = step->GetPreStepPoint();
  auto* mat = pre->GetMaterial();
  auto* pv  = pre->GetPhysicalVolume();
  if (!mat || !pv) return;

  if (track->GetTrackID() == 1 && track->GetCurrentStepNumber() == 1) {
    const auto& v = track->GetVertexPosition();
    fEventAction->SetScanStart(v.y(), v.z());
  }

  const G4double len = step->GetStepLength();
  fEventAction->AddBudget(pv->GetLogicalVolume()->GetName(),
                          len,
                          len / mat->GetRadlen(),
                          len / mat->GetNuclearInterLength(),
                          len * mat->GetDensity());

  if (fStopAtLXe) {
    auto* post = step->GetPostStepPoint()->GetPhysicalVolume();
    if (post && post->GetLogicalVolume()->GetName() == "LXe")
      track->SetTrackStatus(fStopAndKill);
  }
}
