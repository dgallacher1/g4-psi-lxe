#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <vector>

class G4LogicalVolume;

// Geometry (all numbers live in GeometryParameters.hh; frame: z = vessel axis,
// z = 0 at beam-window centre, beam along +x):
//
//   World
//   |- OuterWall / OuterBottom / OuterTop      (steel; wall has the port bore subtracted)
//   |- OuterStub, OuterFlange                  (stub cut by the vessel's outer cylinder)
//   |- AlWindow                                (thin foil in the outer flange bore)
//   `- JacketVac                               (cavity + port bore up to the Al window)
//      |- Sleeve                               (free-standing tube, optional)
//      |- InnerWall / InnerBottom / InnerTop
//      |- InnerStub, InnerFlange, TiWindow
//      `- LXe                                  (cavity + bore up to the Ti window)  <- XenonSD
//         |- SiPM1, SiPM2                      (thin LXe slabs)                     <- PlaneSD
//         `- AlphaSource                       (GDML import)
//
// Every solid is built directly in WORLD coordinates (G4DisplacedSolid), so every
// placement is the identity and the numbers in the header are the numbers in the world.
class DetectorConstruction : public G4VUserDetectorConstruction
{
 public:
  DetectorConstruction() = default;
  ~DetectorConstruction() override = default;

  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;

 private:
  void DefineMaterials();

  G4LogicalVolume* fXenonLV = nullptr;
  std::vector<G4LogicalVolume*> fPlaneLVs;

  G4String fSourceGDML   = "gdml/alpha_source.gdml";
  G4String fSourceLVName = "AlphaSourceLV";
};

#endif
