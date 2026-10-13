#ifndef GeometryBuilders_hh
#define GeometryBuilders_hh 1

// Builds every solid of the cryostat from GeometryParameters.hh.
// Header-only on purpose: DetectorConstruction uses it to build the real
// geometry and tests/geometry_probe.cc uses the very same code to check it.
//
// All solids are expressed in the WORLD frame (G4DisplacedSolid), so booleans
// need no transforms and every placement is the identity.

#include "GeometryParameters.hh"

#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4Transform3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <cmath>

namespace GeoBuild {

// Cylinder/annulus coaxial with z spanning [zMin, zMax].
inline G4VSolid* WorldCyl(const G4String& name, G4double rIn, G4double rOut,
                          G4double zMin, G4double zMax)
{
  auto* t = new G4Tubs(name + "_tubs", rIn, rOut, 0.5 * (zMax - zMin), 0., twopi);
  return new G4DisplacedSolid(name, t, G4Translate3D(0., 0., 0.5 * (zMax + zMin)));
}

// Tube/disc along the beam port spanning s in [s0, s1], s = -x = distance from
// the vessel axis, outward.  rotateY(-90 deg) maps local (x,y,z) -> (-z, y, x),
// so local +z (the tube axis) points along world -x.
inline G4VSolid* PortTube(const G4String& name, G4double rIn, G4double rOut,
                          G4double s0, G4double s1)
{
  auto* t = new G4Tubs(name + "_tubs", rIn, rOut, 0.5 * (s1 - s0), 0., twopi);
  G4RotationMatrix rot;
  rot.rotateY(-90. * deg);
  return new G4DisplacedSolid(name, t, G4Transform3D(rot, G4ThreeVector(-0.5 * (s0 + s1), 0., 0.)));
}

// Start of a stub: its end disc lies wholly inside the vessel's outer
// cylinder, so the cut leaves a clean saddle-shaped weld surface.
inline G4double StubStart(G4double vesselOuterR, G4double stubOR)
{
  return std::sqrt(vesselOuterR * vesselOuterR - stubOR * stubOR) - 1. * mm;
}

struct Solids {
  // outer vessel (placed in World)
  G4VSolid *outerWall, *outerBottom, *outerTop, *outerStub, *outerFlange, *alWindow;
  // insulating vacuum (placed in World); everything below is its child
  G4VSolid *jacket, *sleeve;
  // inner vessel
  G4VSolid *innerWall, *innerBottom, *innerTop, *innerStub, *innerFlange, *tiWindow;
  // liquid xenon (child of jacket)
  G4VSolid *lxe;
  G4double sAlIn, sTiIn;     // inner (vessel-side) faces of the two foils [s]
};

inline Solids Build()
{
  using namespace Geo;
  Solids S{};
  const G4double RO = outerR(), RI = innerR();
  const auto& O = outer;  const auto& I = inner;
  const auto& OP = outerPort;  const auto& IP = innerPort;

  // ---- outer vessel ---------------------------------------------------------
  // Wall = tube minus the port bore (hole radius == stub bore).
  S.outerWall = new G4SubtractionSolid(
      "OuterWall", WorldCyl("OuterWall_cyl", O.cavityR, RO, O.cavityZmin, O.cavityZmax),
      PortTube("OuterWall_hole", 0., OP.boreR, 0., OP.flangeS0));
  S.outerBottom = WorldCyl("OuterBottom", 0., RO, O.cavityZmin - O.bottomT, O.cavityZmin);
  S.outerTop    = WorldCyl("OuterTop",    0., RO, O.cavityZmax, O.cavityZmax + O.topT);
  // Stub = horizontal tube whose root is cut by the vertical vessel's outer surface.
  S.outerStub = new G4SubtractionSolid(
      "OuterStub", PortTube("OuterStub_raw", OP.boreR, OP.stubOR, StubStart(RO, OP.stubOR), OP.flangeS0),
      WorldCyl("OuterStub_cut", 0., RO, -2. * OP.stubOR, 2. * OP.stubOR));
  S.outerFlange = PortTube("OuterFlange", OP.boreR, OP.flangeOR, OP.flangeS0, OP.flangeS1);
  S.sAlIn = alWindow.outerFaceS - alWindow.thickness;
  S.alWindow = PortTube("AlWindow", 0., OP.boreR, S.sAlIn, alWindow.outerFaceS);

  // ---- insulating vacuum: outer cavity + port bore up to the Al window ------
  S.jacket = new G4UnionSolid(
      "JacketVac", WorldCyl("JacketVac_cyl", 0., O.cavityR, O.cavityZmin, O.cavityZmax),
      PortTube("JacketVac_bore", 0., OP.boreR, 0., S.sAlIn));
  S.sleeve = sleeve.use ? PortTube("Sleeve", sleeve.IR, sleeve.OR, sleeve.s0, sleeve.s1) : nullptr;

  // ---- inner vessel ---------------------------------------------------------
  S.innerWall = new G4SubtractionSolid(
      "InnerWall", WorldCyl("InnerWall_cyl", I.cavityR, RI, I.cavityZmin, I.cavityZmax),
      PortTube("InnerWall_hole", 0., IP.boreR, 0., IP.flangeS0));
  S.innerBottom = WorldCyl("InnerBottom", 0., RI, I.cavityZmin - I.bottomT, I.cavityZmin);
  S.innerTop    = WorldCyl("InnerTop",    0., RI, I.cavityZmax, I.cavityZmax + I.topT);
  S.innerStub = new G4SubtractionSolid(
      "InnerStub", PortTube("InnerStub_raw", IP.boreR, IP.stubOR, StubStart(RI, IP.stubOR), IP.flangeS0),
      WorldCyl("InnerStub_cut", 0., RI, -2. * IP.stubOR, 2. * IP.stubOR));
  S.innerFlange = PortTube("InnerFlange", IP.boreR, IP.flangeOR, IP.flangeS0, IP.flangeS1);
  S.sTiIn = tiWindow.outerFaceS - tiWindow.thickness;
  S.tiWindow = PortTube("TiWindow", 0., IP.boreR, S.sTiIn, tiWindow.outerFaceS);

  // ---- liquid xenon: inner cavity + inner port bore up to the Ti foil -------
  S.lxe = new G4UnionSolid(
      "LXe", WorldCyl("LXe_cyl", 0., I.cavityR, I.cavityZmin, I.cavityZmax),
      PortTube("LXe_bore", 0., IP.boreR, 0., S.sTiIn));
  return S;
}

}  // namespace GeoBuild

#endif
