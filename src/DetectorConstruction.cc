#include "DetectorConstruction.hh"
#include "GeometryParameters.hh"
#include "GeometryBuilders.hh"
#include "PlaneSD.hh"
#include "XenonSD.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4GDMLParser.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Transform3D.hh"

#include <fstream>
#include <sstream>

namespace {

G4Material* Mat(const char* name)
{
  auto* m = G4NistManager::Instance()->FindOrBuildMaterial(name);
  if (!m)
    G4Exception("DetectorConstruction", "Geo002", FatalException,
                (G4String("Unknown NIST material: ") + name).c_str());
  return m;
}

}  // namespace

// ------------------------------------------------------------------
void DetectorConstruction::DefineMaterials()
{
  // Liquid xenon is not in the NIST DB with the right condensed-phase state.
  if (!G4Material::GetMaterial("LXe", false)) {
    new G4Material("LXe", /*Z=*/54., /*A=*/131.293 * g / mole,
                   /*rho=*/2.942 * g / cm3, kStateLiquid,
                   /*T=*/165. * kelvin, /*p=*/1.5 * bar);
  }
}

// ------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // ---- refuse to build an incomplete / inconsistent geometry --------------
  const auto problems = Geo::Validate();
  if (!problems.empty()) {
    std::ostringstream os;
    os << "Geometry parameters invalid (" << problems.size() << "):";
    for (const auto& p : problems) os << "\n  - " << p;
    G4Exception("DetectorConstruction::Construct", "Geo001", FatalException, os.str().c_str());
  }
  if (Geo::isDemo) {
    G4cout << "\n"
           << "################################################################\n"
           << "#  DEMO GEOMETRY: placeholder numbers, NOT the CAD model.      #\n"
           << "#  Output is written to xenon_sim_DEMO.root                    #\n"
           << "################################################################\n"
           << G4endl;
  }
  {
    std::ofstream man("geometry_manifest.json");
    man << Geo::ManifestJSON();
  }

  DefineMaterials();
  const G4bool chk = true;

  auto* matWorld  = Mat(Geo::worldMaterial);
  auto* matVac    = Mat("G4_Galactic");
  auto* matOuter  = Mat(Geo::outer.material);
  auto* matInner  = Mat(Geo::inner.material);
  auto* matOPort  = Mat(Geo::outerPort.material);
  auto* matIPort  = Mat(Geo::innerPort.material);
  auto* matAl     = Mat(Geo::alWindow.material);
  auto* matTi     = Mat(Geo::tiWindow.material);
  auto* matSleeve = Mat(Geo::sleeve.material);
  auto* matLXe    = G4Material::GetMaterial("LXe");

  // All solids (world-frame, see GeometryBuilders.hh).
  const GeoBuild::Solids S = GeoBuild::Build();

  auto make = [&](G4VSolid* s, G4Material* m, const G4String& name, G4LogicalVolume* mother,
                  const G4Colour& col) -> G4LogicalVolume* {
    auto* lv = new G4LogicalVolume(s, m, name);
    auto* va = new G4VisAttributes(col);
    va->SetForceSolid(true);
    lv->SetVisAttributes(va);
    new G4PVPlacement(nullptr, G4ThreeVector(), lv, name, mother, false, 0, chk);
    return lv;
  };

  const G4Colour steel(0.6, 0.6, 0.65, 0.5);
  const G4Colour alu(0.85, 0.85, 0.9, 0.8);
  const G4Colour tit(1.0, 0.8, 0.5, 1.0);
  const G4Colour blue(0.2, 0.4, 0.9, 0.6);
  const G4Colour gold(1.0, 0.84, 0.0, 1.0);

  // ---- World ----------------------------------------------------------------
  auto* worldLV = new G4LogicalVolume(
      new G4Box("World", Geo::worldHalfX(), Geo::worldHalfY(), Geo::worldHalfZ()), matWorld, "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV, "World", nullptr, false, 0, chk);

  // ---- outer vessel + Al window (in the World) -------------------------------
  make(S.outerWall,   matOuter, "OuterWall",   worldLV, steel);
  make(S.outerBottom, matOuter, "OuterBottom", worldLV, steel);
  make(S.outerTop,    matOuter, "OuterTop",    worldLV, steel);
  make(S.outerStub,   matOPort, "OuterStub",   worldLV, steel);
  make(S.outerFlange, matOPort, "OuterFlange", worldLV, steel);
  auto* alLV = make(S.alWindow, matAl, "AlWindow", worldLV, alu);

  // ---- insulating vacuum (port bore up to the Al window) ---------------------
  auto* jacketLV = new G4LogicalVolume(S.jacket, matVac, "JacketVac");
  jacketLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  new G4PVPlacement(nullptr, G4ThreeVector(), jacketLV, "JacketVac", worldLV, false, 0, chk);

  if (S.sleeve) make(S.sleeve, matSleeve, "Sleeve", jacketLV, G4Colour(0.1, 0.1, 0.1, 0.8));

  // ---- inner vessel (children of JacketVac) ----------------------------------
  make(S.innerWall,   matInner, "InnerWall",   jacketLV, steel);
  make(S.innerBottom, matInner, "InnerBottom", jacketLV, steel);
  make(S.innerTop,    matInner, "InnerTop",    jacketLV, steel);
  if (S.innerStub) make(S.innerStub, matIPort, "InnerStub", jacketLV, blue);   // none in the current CAD
  make(S.innerFlange, matIPort, "InnerFlange", jacketLV, blue);
  auto* tiLV = make(S.tiWindow, matTi, "TiWindow", jacketLV, tit);

  // ---- liquid xenon (inner cavity + bore up to the Ti foil) -------------------
  fXenonLV = new G4LogicalVolume(S.lxe, matLXe, "LXe");
  auto* lxeVis = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 0.2));
  lxeVis->SetForceSolid(true);
  fXenonLV->SetVisAttributes(lxeVis);
  new G4PVPlacement(nullptr, G4ThreeVector(), fXenonLV, "LXe", jacketLV, false, 0, chk);

  // ---- thin-foil production cuts ---------------------------------------------
  // The default 0.1 mm cut is larger than the foils; give them their own region.
  auto* winRegion = new G4Region("WindowRegion");
  auto* winCuts = new G4ProductionCuts();
  winCuts->SetProductionCut(Geo::windowProductionCut);
  winRegion->SetProductionCuts(winCuts);
  alLV->SetRegion(winRegion); winRegion->AddRootLogicalVolume(alLV);
  tiLV->SetRegion(winRegion); winRegion->AddRootLogicalVolume(tiLV);

  // ---- fluence planes standing in for the SiPMs (silicon slabs inside the LXe) -
  auto* matPlane = Mat(Geo::planeMaterial);
  fPlaneLVs.clear();
  for (G4int i = 0; i < Geo::nPlanes; ++i) {
    const auto& p = Geo::planes[i];
    G4VSolid* shape = p.circular
        ? static_cast<G4VSolid*>(new G4Tubs(p.name, 0., p.halfU, 0.5 * p.thickness, 0., twopi))
        : static_cast<G4VSolid*>(new G4Box(p.name, p.halfU, p.halfV, 0.5 * p.thickness));
    auto* lv = new G4LogicalVolume(shape, matPlane, p.name);
    auto* va = new G4VisAttributes(G4Colour(1.0, 0.3, 0.3, 0.6));
    va->SetForceSolid(true);
    lv->SetVisAttributes(va);
    new G4PVPlacement(G4Transform3D(Geo::PlaneRotation(p), G4ThreeVector(p.cx, p.cy, p.cz)),
                      lv, p.name, fXenonLV, false, /*copyNo=*/i, chk);
    fPlaneLVs.push_back(lv);
  }

  // ---- imported alpha source (GDML), submerged in the LXe ---------------------
  std::ifstream probe(fSourceGDML);
  if (probe.good()) {
    probe.close();
    G4GDMLParser parser;
    parser.SetOverlapCheck(false);
    parser.Read(fSourceGDML, /*validate=*/false);
    if (auto* sourceLV = parser.GetVolume(fSourceLVName)) {
      sourceLV->SetVisAttributes(new G4VisAttributes(gold));
    if (auto* activeLV = parser.GetVolume("ActiveLayerLV"))
      activeLV->SetVisAttributes(new G4VisAttributes(gold));
    if (auto* sourceLV = parser.GetVolume(fSourceLVName)) {
      new G4PVPlacement(nullptr, G4ThreeVector(Geo::sourceX, Geo::sourceY, Geo::sourceZ),
                        sourceLV, "AlphaSource", fXenonLV, false, 0, chk);
    } else {
      G4cerr << "[Detector] WARNING: LV '" << fSourceLVName << "' not in " << fSourceGDML
             << "; running without source." << G4endl;
    }
  }
  } else {
    G4cerr << "[Detector] WARNING: " << fSourceGDML << " not found; running without source." << G4endl;
  }

  G4cout << "[Detector] beam/scan start plane suggestion: x = " << Geo::suggestedBeamStartX() / mm
         << " mm (world half-x = " << Geo::worldHalfX() / mm << " mm)" << G4endl;

  return worldPV;
}

// ------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField()
{
  auto* sdm = G4SDManager::GetSDMpointer();

  // Scoring volume: the inner xenon.
  auto* xenonSD = new XenonSD("XenonSD");
  sdm->AddNewDetector(xenonSD);
  SetSensitiveDetector(fXenonLV, xenonSD);

  // Fluence planes: one SD instance shared by both slabs (plane id = copy number).
  auto* planeSD = new PlaneSD("PlaneSD");
  sdm->AddNewDetector(planeSD);
  for (auto* lv : fPlaneLVs) SetSensitiveDetector(lv, planeSD);
}
