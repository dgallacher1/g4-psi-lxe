#include "EventAction.hh"
#include "XenonHit.hh"
#include "PlaneHit.hh"
#include "OutputOptions.hh"
#include "OutputSchema.hh"

#include "G4Event.hh"
#include "G4PrimaryVertex.hh"
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

// ---- one row per LXe hit, one function per ntuple ---------------------------------------------

void FillStepRow(G4AnalysisManager* man, G4int eid, const XenonHit& h)
{
  using N = Out::Steps;
  Out::Row<N> r(man);
  r.set(N::eventID, eid);
  r.set(N::trackID, h.trackID);
  r.set(N::parentID, h.parentID);
  r.set(N::pdg, h.pdg);
  r.set(N::preX_mm, h.prePos.x() / mm);
  r.set(N::preY_mm, h.prePos.y() / mm);
  r.set(N::preZ_mm, h.prePos.z() / mm);
  r.set(N::postX_mm, h.postPos.x() / mm);
  r.set(N::postY_mm, h.postPos.y() / mm);
  r.set(N::postZ_mm, h.postPos.z() / mm);
  r.set(N::edep_keV, h.edep / keV);
  r.set(N::kinE_keV, h.kinE / keV);
  r.set(N::stepLen_mm, h.stepLen / mm);
  r.set(N::time_ns, h.gtime / ns);
  r.set(N::niel_keV, h.niel / keV);
  r.add();
}

void FillChromaRow(G4AnalysisManager* man, G4int eid, G4int row, const XenonHit& h)
{
  using N = Out::Chroma;
  Out::Row<N> r(man);
  r.set(N::event, eid);
  r.set(N::hits, row);
  r.set(N::xpos, h.postPos.x() / mm);                 // post-step point
  r.set(N::ypos, h.postPos.y() / mm);
  r.set(N::zpos, h.postPos.z() / mm);
  r.set(N::time, h.postTime / ns);
  r.set(N::edep, h.edep / keV);
  r.set(N::energy, h.kinE / keV);                     // pre-step kinetic energy
  r.set(N::x_momentum, h.preMom.x() / keV);           // pre-step momentum
  r.set(N::y_momentum, h.preMom.y() / keV);
  r.set(N::z_momentum, h.preMom.z() / keV);
  r.set(N::particle_name, ParticleClass(h.pdg));
  r.set(N::PDG, h.pdg);
  r.set(N::trackID, h.trackID);
  r.add();
}

void FillFuseRow(G4AnalysisManager* man, G4int eid, const XenonHit& h, const G4ThreeVector& vtx)
{
  using N = Out::Fuse;
  Out::Row<N> r(man);
  const TrackLabel* L = h.label;                      // null only if TrackingAction was not registered
  const G4double x = h.postPos.x() / mm, y = h.postPos.y() / mm, z = h.postPos.z() / mm;
  const G4double xv = vtx.x() / mm, yv = vtx.y() / mm, zv = vtx.z() / mm;
  r.set(N::eventid, eid);
  r.set(N::trackid, h.trackID);
  r.set(N::parentid, h.parentID);
  r.set(N::type, L ? L->type : std::string("unknown"));
  r.set(N::parenttype, L ? L->parenttype : std::string("none"));
  r.set(N::creaproc, L ? L->creaproc : std::string("none"));
  r.set(N::edproc, h.edproc);
  r.set(N::ed, h.edep / keV);
  r.set(N::t, h.postTime / ns);
  r.set(N::xp, x);  r.set(N::yp, y);  r.set(N::zp, z);
  r.set(N::xp_pri, xv);  r.set(N::yp_pri, yv);  r.set(N::zp_pri, zv);
  r.set(N::x, x);  r.set(N::y, y);  r.set(N::z, z);                      // duplicates, see OutputSchema.hh
  r.set(N::x_pri, xv);  r.set(N::y_pri, yv);  r.set(N::z_pri, zv);
  r.add();
}

}  // namespace

// -----------------------------------------------------------------------------------------------

void EventAction::BeginOfEventAction(const G4Event*)
{
  fPath.clear();
  gTrackLabels.clear();               // per-event track labels (filled by TrackingAction)
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
  const G4int eid = evt->GetEventID();

  G4ThreeVector primaryVertex;
  if (auto* pv = evt->GetPrimaryVertex()) primaryVertex = G4ThreeVector(pv->GetX0(), pv->GetY0(), pv->GetZ0());

  Totals tot;
  if (auto* hce = evt->GetHCofThisEvent()) {
    FillLXeHits(hce, eid, primaryVertex, tot);
    FillPlaneHits(hce, eid, tot);
  }
  FillPrimaryPath(eid);
  FillEventSummary(eid, tot);
  FillBudget(eid);
}

// One pass over the LXe hits: steps (all hits), tree_chroma (edep > 0), tree_fuse (edep > 0 and the
// gamma interaction steps, which the lineage clustering needs as parents).
void EventAction::FillLXeHits(G4HCofThisEvent* hce, G4int eid, const G4ThreeVector& primaryVertex, Totals& tot)
{
  auto* man = G4AnalysisManager::Instance();
  if (fXenonHCID < 0) fXenonHCID = G4SDManager::GetSDMpointer()->GetCollectionID("XenonHits");
  auto* hc = static_cast<XenonHitsCollection*>(hce->GetHC(fXenonHCID));
  if (!hc) return;

  tot.nSteps = hc->entries();
  G4int chromaRow = 0;                                  // row index within this event in tree_chroma
  for (G4int i = 0; i < tot.nSteps; ++i) {
    const XenonHit& h = *(*hc)[i];
    tot.lxeEdep += h.edep;

    if (gWriteChroma && h.edep > 0.) FillChromaRow(man, eid, chromaRow++, h);
    if (gWriteFuse && (h.edep > 0. || (h.pdg == 22 && h.edproc != "Transportation")))
      FillFuseRow(man, eid, h, primaryVertex);
    if (gWriteSteps) FillStepRow(man, eid, h);
  }
}

// Charged-particle steps in the fluence planes.
void EventAction::FillPlaneHits(G4HCofThisEvent* hce, G4int eid, Totals& tot)
{
  using N = Out::Fluence;
  auto* man = G4AnalysisManager::Instance();
  if (fPlaneHCID < 0) fPlaneHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PlaneHits");
  auto* hc = static_cast<PlaneHitsCollection*>(hce->GetHC(fPlaneHCID));
  if (!hc) return;

  Out::Row<N> r(man);
  const G4int nPlane = hc->entries();
  for (G4int i = 0; i < nPlane; ++i) {
    const PlaneHit& h = *(*hc)[i];
    r.set(N::eventID, eid);
    r.set(N::planeID, h.planeID);
    r.set(N::trackID, h.trackID);
    r.set(N::parentID, h.parentID);
    r.set(N::pdg, h.pdg);
    r.set(N::entering, h.entering);
    r.set(N::u_mm, h.u / mm);
    r.set(N::v_mm, h.v / mm);
    r.set(N::kinE_keV, h.kinE / keV);
    r.set(N::cosTheta, h.cosTheta);
    r.set(N::stepLen_mm, h.stepLen / mm);
    r.set(N::edep_keV, h.edep / keV);
    r.set(N::time_ns, h.gtime / ns);
    r.set(N::vtxX_mm, h.vx / mm);
    r.set(N::vtxY_mm, h.vy / mm);
    r.set(N::vtxZ_mm, h.vz / mm);
    r.set(N::creator, h.creator);
    r.add();
    tot.planeEdep += h.edep;
  }
}

// The beam particle's start point and every step end, in every volume (filled by SteppingAction).
void EventAction::FillPrimaryPath(G4int eid) const
{
  if (!gWritePath) return;
  using N = Out::Path;
  Out::Row<N> r(G4AnalysisManager::Instance());
  G4int k = 0;
  for (const auto& q : fPath) {
    r.set(N::eventID, eid);
    r.set(N::step, k++);
    r.set(N::x_mm, q.x / mm);
    r.set(N::y_mm, q.y / mm);
    r.set(N::z_mm, q.z / mm);
    r.set(N::kinE_keV, q.kinE / keV);
    r.set(N::time_ns, q.t / ns);
    r.set(N::volume, q.vol);
    r.add();
  }
}

// edep_keV is the energy deposited in the xenon; edepPlanes_keV the energy deposited in the silicon
// planes (kept separate: it is not xenon signal).
void EventAction::FillEventSummary(G4int eid, const Totals& tot) const
{
  using N = Out::Events;
  Out::Row<N> r(G4AnalysisManager::Instance());
  r.set(N::eventID, eid);
  r.set(N::edep_keV, tot.lxeEdep / keV);
  r.set(N::nSteps, tot.nSteps);
  r.set(N::edepPlanes_keV, tot.planeEdep / keV);
  r.add();
}

// Material budget per volume (geantino scans only: fBudget is empty otherwise).
void EventAction::FillBudget(G4int eid) const
{
  using N = Out::Budget;
  Out::Row<N> r(G4AnalysisManager::Instance());
  for (const auto& kv : fBudget) {
    r.set(N::eventID, eid);
    r.set(N::volume, kv.first);
    r.set(N::y0_mm, fY0 / mm);
    r.set(N::z0_mm, fZ0 / mm);
    r.set(N::length_mm, kv.second.len / mm);
    r.set(N::X0, kv.second.x0);
    r.set(N::lambdaI, kv.second.lambda);
    r.set(N::massThick_gcm2, kv.second.massThick / (g / cm2));
    r.add();
  }
}
