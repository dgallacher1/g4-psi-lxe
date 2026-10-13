#include "EventAction.hh"
#include "XenonHit.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

void EventAction::EndOfEventAction(const G4Event* evt)
{
  auto* hce = evt->GetHCofThisEvent();
  if (!hce) return;

  if (fHCID < 0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID("XenonHits");

  auto* hc = static_cast<XenonHitsCollection*>(hce->GetHC(fHCID));
  if (!hc) return;

  auto* man = G4AnalysisManager::Instance();
  const G4int eid = evt->GetEventID();
  const G4int n   = hc->entries();
  G4double totalEdep = 0.;

  for (G4int i = 0; i < n; ++i) {
    auto* h = (*hc)[i];

    man->FillNtupleIColumn(0, 0, eid);
    man->FillNtupleIColumn(0, 1, h->trackID);
    man->FillNtupleIColumn(0, 2, h->parentID);
    man->FillNtupleIColumn(0, 3, h->pdg);
    man->FillNtupleDColumn(0, 4, h->prePos.x() / mm);
    man->FillNtupleDColumn(0, 5, h->prePos.y() / mm);
    man->FillNtupleDColumn(0, 6, h->prePos.z() / mm);
    man->FillNtupleDColumn(0, 7, h->postPos.x() / mm);
    man->FillNtupleDColumn(0, 8, h->postPos.y() / mm);
    man->FillNtupleDColumn(0, 9, h->postPos.z() / mm);
    man->FillNtupleDColumn(0, 10, h->edep / keV);
    man->FillNtupleDColumn(0, 11, h->kinE / keV);
    man->FillNtupleDColumn(0, 12, h->stepLen / mm);
    man->FillNtupleDColumn(0, 13, h->gtime / ns);
    man->AddNtupleRow(0);

    totalEdep += h->edep;
  }

  man->FillNtupleIColumn(1, 0, eid);
  man->FillNtupleDColumn(1, 1, totalEdep / keV);
  man->FillNtupleIColumn(1, 2, n);
  man->AddNtupleRow(1);
}
