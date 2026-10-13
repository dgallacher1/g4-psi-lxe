#include "PlaneSD.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

PlaneSD::PlaneSD(const G4String& name) : G4VSensitiveDetector(name)
{
  collectionName.insert("PlaneHits");
  for (G4int i = 0; i < Geo::nPlanes; ++i) fRot[i] = Geo::PlaneRotation(Geo::planes[i]);
}

void PlaneSD::Initialize(G4HCofThisEvent* hce)
{
  fHC = new PlaneHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHCID < 0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(fHCID, fHC);
}

G4bool PlaneSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto* track = step->GetTrack();
  if (track->GetDefinition()->GetPDGCharge() == 0.) return false;   // charged-particle fluence only

  auto* pre = step->GetPreStepPoint();
  const G4int id = pre->GetTouchable()->GetCopyNumber();
  if (id < 0 || id >= Geo::nPlanes) return false;

  const auto& spec = Geo::planes[id];
  const G4RotationMatrix& R = fRot[id];
  const G4ThreeVector centre(spec.cx, spec.cy, spec.cz);
  const G4ThreeVector local  = R.inverse() * (pre->GetPosition() - centre);
  const G4ThreeVector normal = R * G4ThreeVector(0., 0., 1.);

  auto* hit = new PlaneHit();
  hit->planeID  = id;
  hit->trackID  = track->GetTrackID();
  hit->parentID = track->GetParentID();
  hit->pdg      = track->GetDefinition()->GetPDGEncoding();
  hit->entering = (pre->GetStepStatus() == fGeomBoundary) ? 1 : 0;
  hit->u        = local.x();
  hit->v        = local.y();
  hit->kinE     = pre->GetKineticEnergy();
  hit->cosTheta = pre->GetMomentumDirection().dot(normal);
  hit->stepLen  = step->GetStepLength();
  hit->edep     = step->GetTotalEnergyDeposit();
  hit->gtime    = pre->GetGlobalTime();
  hit->vx       = track->GetVertexPosition().x();
  hit->vy       = track->GetVertexPosition().y();
  hit->vz       = track->GetVertexPosition().z();
  hit->creator  = track->GetCreatorProcess() ? track->GetCreatorProcess()->GetProcessName() : "primary";

  fHC->insert(hit);
  return true;
}
