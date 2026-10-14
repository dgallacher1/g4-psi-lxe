#include "TrackingAction.hh"
#include "TrackLabels.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ParticleDefinition.hh"

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  TrackLabel label;
  label.type = track->GetDefinition()->GetParticleName();

  label.parenttype = "none";
  if (track->GetParentID() > 0) {
    // the parent is always tracked (and labelled) before its secondaries
    auto it = gTrackLabels.find(track->GetParentID());
    if (it != gTrackLabels.end()) label.parenttype = it->second.type;
  }

  const G4VProcess* creator = track->GetCreatorProcess();
  label.creaproc = creator ? creator->GetProcessName() : "none";
  // the lineage code (Clustering.py) knows the decay process under this name
  if (label.creaproc.rfind("RadioactiveDecay", 0) == 0) label.creaproc = "RadioactiveDecayBase";

  gTrackLabels[track->GetTrackID()] = std::move(label);
}
