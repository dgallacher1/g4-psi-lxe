#include "EventAction.hh"
#include "XenonHit.hh"
#include "PlaneHit.hh"
#include "OutputOptions.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

namespace {
// numeric class written to tree_chroma/particle_name: 1 gamma, 2 neutron, 3 e-, 4 e+, 5 other
G4double ParticleClass(G4int pdg)
{
  switch (pdg) {
    case 22:   return 1.;
    case 2112: return 2.;
    case 11:   return 3.;
    case -11:  return 4.;
    default:   return 5.;
  }
}
}  // namespace

void EventAction::BeginOfEventAction(const G4Event*)
{
  fBudget.clear();
  fY0 = fZ0 = 0.;
}

void EventAction::AddBudget(const G4String& volume, G4double len, G4double x0,
                            G4double lambda, G4double massThick)
{
  auto& b = fBudget[volume];
  b.len += len; b.x0 += x0; b.lambda += lambda; b.massThick += massThick;
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  auto* man = G4AnalysisManager::Instance();
  const G4int eid = evt->GetEventID();

  G4double totalEdep = 0., planeEdep = 0.;
  G4int nSteps = 0;
  G4int chromaRow = 0;                // row index within this event in tree_chroma

  if (auto* hce = evt->GetHCofThisEvent()) {
    auto* sdm = G4SDManager::GetSDMpointer();
    if (fXenonHCID < 0) fXenonHCID = sdm->GetCollectionID("XenonHits");
    if (fPlaneHCID < 0) fPlaneHCID = sdm->GetCollectionID("PlaneHits");

    // ---- ntuple 0: every step in the LXe ------------------------------------
    if (auto* hc = static_cast<XenonHitsCollection*>(hce->GetHC(fXenonHCID))) {
      nSteps = hc->entries();
      for (G4int i = 0; i < nSteps; ++i) {
        auto* h = (*hc)[i];
        totalEdep += h->edep;

        if (gWriteChroma && h->edep > 0.) {          // tree_chroma: energy-depositing steps only
          man->FillNtupleDColumn(4, 0, eid);
          man->FillNtupleDColumn(4, 1, chromaRow++);
          man->FillNtupleDColumn(4, 2, h->postPos.x() / mm);      // post-step point
          man->FillNtupleDColumn(4, 3, h->postPos.y() / mm);
          man->FillNtupleDColumn(4, 4, h->postPos.z() / mm);
          man->FillNtupleDColumn(4, 5, h->postTime / ns);
          man->FillNtupleDColumn(4, 6, h->edep / keV);
          man->FillNtupleDColumn(4, 7, h->kinE / keV);            // pre-step kinetic energy
          man->FillNtupleDColumn(4, 8, h->preMom.x() / keV);      // pre-step momentum
          man->FillNtupleDColumn(4, 9, h->preMom.y() / keV);
          man->FillNtupleDColumn(4, 10, h->preMom.z() / keV);
          man->FillNtupleDColumn(4, 11, ParticleClass(h->pdg));
          man->FillNtupleDColumn(4, 12, h->pdg);
          man->FillNtupleDColumn(4, 13, h->trackID);
          man->AddNtupleRow(4);
        }

        if (!gWriteSteps) continue;
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
        man->FillNtupleDColumn(0, 14, h->niel / keV);
        man->AddNtupleRow(0);
      }
    }

    // ---- ntuple 2: charged-particle steps in the fluence planes -------------
    if (auto* hc = static_cast<PlaneHitsCollection*>(hce->GetHC(fPlaneHCID))) {
      const G4int n = hc->entries();
      for (G4int i = 0; i < n; ++i) {
        auto* h = (*hc)[i];
        man->FillNtupleIColumn(2, 0, eid);
        man->FillNtupleIColumn(2, 1, h->planeID);
        man->FillNtupleIColumn(2, 2, h->trackID);
        man->FillNtupleIColumn(2, 3, h->parentID);
        man->FillNtupleIColumn(2, 4, h->pdg);
        man->FillNtupleIColumn(2, 5, h->entering);
        man->FillNtupleDColumn(2, 6, h->u / mm);
        man->FillNtupleDColumn(2, 7, h->v / mm);
        man->FillNtupleDColumn(2, 8, h->kinE / keV);
        man->FillNtupleDColumn(2, 9, h->cosTheta);
        man->FillNtupleDColumn(2, 10, h->stepLen / mm);
        man->FillNtupleDColumn(2, 11, h->edep / keV);
        man->FillNtupleDColumn(2, 12, h->gtime / ns);
        man->FillNtupleDColumn(2, 13, h->vx / mm);
        man->FillNtupleDColumn(2, 14, h->vy / mm);
        man->FillNtupleDColumn(2, 15, h->vz / mm);
        man->FillNtupleSColumn(2, 16, h->creator);
        man->AddNtupleRow(2);
        planeEdep += h->edep;
      }
    }
  }

  // ---- ntuple 1: per-event summary -------------------------------------------
  // edep_keV is the energy deposited in the xenon; edepPlanes_keV is the energy
  // deposited in the silicon planes (kept separate: it is not xenon signal).
  man->FillNtupleIColumn(1, 0, eid);
  man->FillNtupleDColumn(1, 1, totalEdep / keV);
  man->FillNtupleIColumn(1, 2, nSteps);
  man->FillNtupleDColumn(1, 3, planeEdep / keV);
  man->AddNtupleRow(1);

  // ---- ntuple 3: material budget per volume (geantino scans) -----------------
  for (const auto& kv : fBudget) {
    man->FillNtupleIColumn(3, 0, eid);
    man->FillNtupleSColumn(3, 1, kv.first);
    man->FillNtupleDColumn(3, 2, fY0 / mm);
    man->FillNtupleDColumn(3, 3, fZ0 / mm);
    man->FillNtupleDColumn(3, 4, kv.second.len / mm);
    man->FillNtupleDColumn(3, 5, kv.second.x0);
    man->FillNtupleDColumn(3, 6, kv.second.lambda);
    man->FillNtupleDColumn(3, 7, kv.second.massThick / (g / cm2));
    man->AddNtupleRow(3);
  }
}
