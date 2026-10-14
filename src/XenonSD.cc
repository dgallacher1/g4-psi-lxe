#include "XenonSD.hh"

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4GenericMessenger.hh"
#include "G4VProcess.hh"
#include "G4GammaGeneralProcess.hh"
#include "G4Track.hh"

XenonSD::XenonSD(const G4String& name) : G4VSensitiveDetector(name)
{
  collectionName.insert("XenonHits");
  fMessenger = new G4GenericMessenger(this, "/lxe/", "LXe sensitive-detector control");
  fMessenger->DeclareProperty("onlyEdep", fOnlyEdep,
      "true: store only steps with edep > 0 (default false: every step, incl. neutral fly-throughs)");
}

XenonSD::~XenonSD() { delete fMessenger; }

void XenonSD::Initialize(G4HCofThisEvent* hce)
{
  fHC = new XenonHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHCID < 0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(fHCID, fHC);
}

G4bool XenonSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  // Called for every step whose PRE-step point is in the LXe. We record all
  // steps (not just edep>0) so the trajectory polyline is complete, including
  // neutral fly-throughs. To decimate, gate on edep or on charge here.
  auto* pre   = step->GetPreStepPoint();
  auto* track = step->GetTrack();

  // Process that ended the step. With the gamma general process active the sub-process
  // ("compt", "phot", "conv", "Rayl") is what identifies the interaction.
  std::string edproc = "none";
  if (const G4VProcess* proc = step->GetPostStepPoint()->GetProcessDefinedStep()) {
    edproc = proc->GetProcessName();
    if (track->GetDefinition()->GetPDGEncoding() == 22) {
      if (auto* gg = dynamic_cast<const G4GammaGeneralProcess*>(proc)) edproc = gg->GetSubProcessName();
    }
  }

  // /lxe/onlyEdep keeps energy deposits and a gamma's discrete interactions (the lineage
  // clustering needs the gamma rows as parents of the electrons they release).
  const G4bool gammaInteraction = (track->GetDefinition()->GetPDGEncoding() == 22) && (edproc != "Transportation");
  if (fOnlyEdep && step->GetTotalEnergyDeposit() <= 0. && !gammaInteraction) return false;

  auto* hit = new XenonHit();
  hit->trackID  = track->GetTrackID();
  hit->parentID = track->GetParentID();
  hit->pdg      = track->GetDefinition()->GetPDGEncoding();
  hit->edep     = step->GetTotalEnergyDeposit();
  hit->niel     = step->GetNonIonizingEnergyDeposit();
  hit->kinE     = pre->GetKineticEnergy();
  hit->stepLen  = step->GetStepLength();
  hit->gtime    = pre->GetGlobalTime();
  hit->postTime = step->GetPostStepPoint()->GetGlobalTime();
  hit->preMom   = pre->GetMomentum();
  hit->edproc   = std::move(edproc);
  auto lit = gTrackLabels.find(track->GetTrackID());
  hit->label    = (lit != gTrackLabels.end()) ? &lit->second : nullptr;
  hit->prePos   = pre->GetPosition();
  hit->postPos  = step->GetPostStepPoint()->GetPosition();

  fHC->insert(hit);
  return true;
}
