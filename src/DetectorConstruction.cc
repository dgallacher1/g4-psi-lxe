#include "DetectorConstruction.hh"
#include "XenonSD.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GDMLParser.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <fstream>

// ------------------------------------------------------------------
void DetectorConstruction::DefineMaterials()
{
  auto* nist = G4NistManager::Instance();
  nist->FindOrBuildMaterial("G4_Galactic");          // world + vacuum jacket
  nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");   // vessels

  // Liquid xenon is NOT in the NIST DB with the right condensed-phase
  // density/state, so build it explicitly. ~2.94 g/cm^3 near 165 K.
  if (!G4Material::GetMaterial("LXe", false)) {
    new G4Material("LXe", /*Z=*/54., /*A=*/131.293 * g / mole,
                   /*rho=*/2.942 * g / cm3, kStateLiquid,
                   /*T=*/165. * kelvin, /*p=*/1.5 * bar);
  }
}

// ------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();

  auto* vac   = G4Material::GetMaterial("G4_Galactic");
  auto* steel = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto* lxe   = G4Material::GetMaterial("LXe");

  const G4bool checkOverlaps = true;

  // ---- radial/axial build parameters (half-lengths where noted) ----
  const G4double rLXe = 150. * mm, hzLXe = 150. * mm;   // active LXe
  const G4double tInner = 3. * mm;                      // inner vessel wall
  const G4double tGap   = 30. * mm;                     // vacuum jacket
  const G4double tOuter = 5. * mm;                      // outer vessel wall

  const G4double rInner = rLXe + tInner,  hzInner = hzLXe + tInner;
  const G4double rGap   = rInner + tGap,  hzGap   = hzInner + tGap;
  const G4double rOuter = rGap + tOuter,  hzOuter = hzGap + tOuter;
  const G4double worldHalf = rOuter + 100. * mm;

  // ---- World ----
  auto* worldSolid = new G4Box("World", worldHalf, worldHalf, worldHalf);
  auto* worldLV    = new G4LogicalVolume(worldSolid, vac, "World");
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "World",
                                    nullptr, false, 0, checkOverlaps);

  // Nested solid cylinders: each daughter "carves" the parent, so the steel
  // vessels are automatically closed (walls + end caps) without booleans.
  auto* outerSolid = new G4Tubs("OuterVessel", 0, rOuter, hzOuter, 0, twopi);
  auto* outerLV = new G4LogicalVolume(outerSolid, steel, "OuterVessel");
  new G4PVPlacement(nullptr, {}, outerLV, "OuterVessel", worldLV, false, 0, checkOverlaps);

  auto* gapSolid = new G4Tubs("VacuumJacket", 0, rGap, hzGap, 0, twopi);
  auto* gapLV = new G4LogicalVolume(gapSolid, vac, "VacuumJacket");
  new G4PVPlacement(nullptr, {}, gapLV, "VacuumJacket", outerLV, false, 0, checkOverlaps);

  auto* innerSolid = new G4Tubs("InnerVessel", 0, rInner, hzInner, 0, twopi);
  auto* innerLV = new G4LogicalVolume(innerSolid, steel, "InnerVessel");
  new G4PVPlacement(nullptr, {}, innerLV, "InnerVessel", gapLV, false, 0, checkOverlaps);

  auto* lxeSolid = new G4Tubs("LXe", 0, rLXe, hzLXe, 0, twopi);
  fXenonLV = new G4LogicalVolume(lxeSolid, lxe, "LXe");
  auto* lxeVis = new G4VisAttributes(G4Colour(0.2, 0.6, 1.0, 0.25));
  lxeVis->SetForceSolid(true);
  fXenonLV->SetVisAttributes(lxeVis);
  new G4PVPlacement(nullptr, {}, fXenonLV, "LXe", innerLV, false, 0, checkOverlaps);

  // ---- Imported alpha source (GDML), submerged at LXe centre ----
  std::ifstream probe(fSourceGDML);
  if (probe.good()) {
    probe.close();
    G4GDMLParser parser;
    parser.SetOverlapCheck(false);        // tessellated CAD solids: overlap check is costly
    parser.Read(fSourceGDML, /*validate=*/false);  // CAD->GDML often not schema-perfect
    auto* sourceLV = parser.GetVolume(fSourceLVName);
    if (sourceLV) {
      new G4PVPlacement(nullptr, G4ThreeVector(0, 0, 0), sourceLV,
                        "AlphaSource", fXenonLV, false, 0, checkOverlaps);
      G4cout << "[Detector] Placed GDML source '" << fSourceLVName
             << "' inside LXe." << G4endl;
    } else {
      G4cerr << "[Detector] WARNING: LV '" << fSourceLVName
             << "' not found in " << fSourceGDML << " — running without source."
             << G4endl;
    }
  } else {
    G4cerr << "[Detector] WARNING: " << fSourceGDML
           << " not found — running without the imported source." << G4endl;
  }

  return worldPV;
}

// ------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField()
{
  // One SD on the LXe LV: the scoring volume is exactly the inner xenon.
  auto* sd = new XenonSD("XenonSD");
  G4SDManager::GetSDMpointer()->AddNewDetector(sd);
  SetSensitiveDetector(fXenonLV, sd);
}
