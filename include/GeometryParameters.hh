#ifndef GeometryParameters_hh
#define GeometryParameters_hh 1

// ============================================================================
//  GeometryParameters.hh  --  single source of truth for the cryostat geometry
// ============================================================================
//
//  FRAME (matches the GPU light simulations)
//    +z : long axis of the cylindrical vessels.  z = 0 is the CENTRE OF THE
//         BEAM WINDOW (i.e. the plane containing the port axis).
//    +x : beam direction.  The beam port is on the -x side of the vessels and
//         the beam enters from -x, travelling towards +x.
//    x = y = 0 : the common vessel axis (outer and inner vessels are coaxial).
//    y : completes the right-handed frame (y = z cross x).
//
//  PORT COORDINATE  s
//    Everything on the beam port is located by s = distance from the vessel
//    axis, measured OUTWARD along the port axis (s = -x).  Larger s = more
//    upstream.  Port parts are therefore listed in the order the beam meets
//    them when you read s from large to small:
//        outer CF flange  ->  Al window  ->  sleeve  ->  inner CF flange / Ti
//        window  ->  LXe.
//
//  HOW TO FILL THIS IN
//    Every number has a GEO_V(<placeholder>) wrapper.  Overwrite the whole
//    GEO_V(...) with the CAD number (plain `12.7*mm`).  Anything you have not
//    replaced is a placeholder demo value:
//      * built with -DGEO_DEMO=ON  : the placeholders are used, the run is
//        loudly flagged and the output file is named xenon_sim_DEMO.root.
//      * built with -DGEO_DEMO=OFF : unreplaced entries are NaN and
//        DetectorConstruction aborts, listing exactly which ones are missing.
//    A CAD-derived geometry can therefore never be silently mixed with
//    placeholders.
//
//  Dimensions are mm (Geant4 internal units).  Material names are NIST names.
// ============================================================================

#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace Geo {

inline constexpr G4double TBD = std::numeric_limits<G4double>::quiet_NaN();

#ifdef GEO_USE_DEMO_VALUES
  inline constexpr G4bool isDemo = true;
  #define GEO_V(demo) (demo)
#else
  inline constexpr G4bool isDemo = false;
  #define GEO_V(demo) (::Geo::TBD)
#endif

// ----------------------------------------------------------------------------
//  Types
// ----------------------------------------------------------------------------
struct Vessel {                       // one cylindrical vessel, coaxial with z
  G4double cavityR;                   // inner radius of the wall
  G4double cavityZmin, cavityZmax;    // inner cavity extent in WORLD z
  G4double wallT;                     // side-wall thickness
  G4double bottomT, topT;             // end-plate thicknesses (added outside the cavity)
  const char* material;
};

struct Port {                         // CF stub + flange, along the port axis
  G4double boreR;                     // through-bore radius (stub ID/2 == flange bore)
  G4double stubOR;                    // stub outer radius
  G4double flangeOR;                  // flange outer radius (OD/2)
  G4double flangeS0, flangeS1;        // flange inner / outer face (s)
  const char* material;
};

struct Window {                       // thin foil filling the port bore
  G4double thickness;
  G4double outerFaceS;                // s of the foil face that points away from the vessel
  const char* material;
};

struct Sleeve {                       // free-standing tube in the vacuum jacket (not attached)
  G4bool   use;
  G4double IR, OR;
  G4double s0, s1;                    // s0 < s1
  const char* material;
};

struct PlaneSpec {                    // thin virtual slab for charged-particle fluence
  const char* name;
  G4double cx, cy, cz;                // centre, world mm
  // Orientation: local frame rotated by rotX, then rotY, then rotZ (degrees,
  // active, about FIXED world axes).  Local z is the plane NORMAL; halfU is the
  // half-extent along local x, halfV along local y.  Common cases:
  //   normal +x : rotY = +90   (U -> -z, V -> +y)
  //   normal -x : rotY = -90   (U -> +z, V -> +y)
  //   normal +z : all zero     (U -> +x, V -> +y)
  //   normal +y : rotX = -90   (U -> +x, V -> -z)
  G4double rotXdeg, rotYdeg, rotZdeg;
  G4double halfU, halfV;
  G4double thickness;                 // keep thin (tens of um) so it does not perturb tracks
};

// ----------------------------------------------------------------------------
//  1. Vessels   (CAD: cut-away section, vertical cylinders)
// ----------------------------------------------------------------------------
inline constexpr Vessel outer = {
  /*cavityR   */ GEO_V(100.*mm),
  /*cavityZmin*/ GEO_V(-230.*mm),
  /*cavityZmax*/ GEO_V(300.*mm),
  /*wallT     */ GEO_V(3.*mm),
  /*bottomT   */ GEO_V(10.*mm),
  /*topT      */ GEO_V(10.*mm),
  "G4_STAINLESS-STEEL"
};

inline constexpr Vessel inner = {
  /*cavityR   */ GEO_V(45.*mm),
  /*cavityZmin*/ GEO_V(-170.*mm),
  /*cavityZmax*/ GEO_V(120.*mm),
  /*wallT     */ GEO_V(2.*mm),
  /*bottomT   */ GEO_V(6.*mm),
  /*topT      */ GEO_V(6.*mm),
  "G4_STAINLESS-STEEL"
};

// ----------------------------------------------------------------------------
//  2. Beam port on the OUTER vessel  (standard DN40CF = 2.75" Conflat)
//     stub is welded to the wall; flange is a simple additive ring.
//     Placeholders below are NOMINAL DN40CF: flange OD 69.85 mm, tube OD
//     38.1 mm (1.5"), ~1.5 mm wall (ID ~35 mm), flange ~15.75 mm thick.
//     Replace with the CAD values.  Knife edge, gasket groove and bolt holes
//     are not modelled (they are away from the beam axis).
// ----------------------------------------------------------------------------
inline constexpr Port outerPort = {
  /*boreR   */ GEO_V(17.5*mm),
  /*stubOR  */ GEO_V(19.05*mm),
  /*flangeOR*/ GEO_V(34.925*mm),
  /*flangeS0*/ GEO_V(110.*mm),
  /*flangeS1*/ GEO_V(125.75*mm),
  "G4_STAINLESS-STEEL"
};

// First (thin) window: aluminium, separates the jacket vacuum from the beamline side.
inline constexpr Window alWindow = {
  /*thickness */ GEO_V(0.1*mm),
  /*outerFaceS*/ GEO_V(125.75*mm),
  "G4_Al"
};

// ----------------------------------------------------------------------------
//  3. Sleeve between the two windows (the dark tube in the CAD image)
// ----------------------------------------------------------------------------
inline constexpr Sleeve sleeve = {
  /*use*/ true,
  /*IR */ GEO_V(15.*mm),
  /*OR */ GEO_V(16.5*mm),
  /*s0 */ GEO_V(70.*mm),
  /*s1 */ GEO_V(120.*mm),
  "G4_Al"
};

// ----------------------------------------------------------------------------
//  4. Port on the INNER vessel: non-standard CF, concentric with the outer
//     port, NOT attached to it (cryo insulation gap).  The Ti foil closes the
//     bore; LXe fills the bore up to the foil.
// ----------------------------------------------------------------------------
inline constexpr Port innerPort = {
  /*boreR   */ GEO_V(17.5*mm),
  /*stubOR  */ GEO_V(20.*mm),
  /*flangeOR*/ GEO_V(34.925*mm),
  /*flangeS0*/ GEO_V(50.*mm),
  /*flangeS1*/ GEO_V(62.7*mm),
  "G4_STAINLESS-STEEL"
};

inline constexpr Window tiWindow = {
  /*thickness */ GEO_V(50.*um),
  /*outerFaceS*/ GEO_V(62.7*mm),
  "G4_Ti"
};

// ----------------------------------------------------------------------------
//  5. Virtual SiPM fluence planes (LXe-filled slabs, children of the LXe)
// ----------------------------------------------------------------------------
inline constexpr G4int nPlanes = 2;

inline constexpr PlaneSpec planes[nPlanes] = {
  // name    cx            cy           cz            rotX         rotY        rotZ        halfU          halfV         thick
  {"SiPM1", GEO_V(-20.*mm), GEO_V(0.*mm), GEO_V(0.*mm),   GEO_V(0.), GEO_V(90.), GEO_V(0.), GEO_V(12.*mm), GEO_V(12.*mm), GEO_V(0.05*mm)},
  {"SiPM2", GEO_V(-10.*mm), GEO_V(0.*mm), GEO_V(-25.*mm), GEO_V(0.), GEO_V(0.),  GEO_V(0.), GEO_V(20.*mm), GEO_V(20.*mm), GEO_V(0.05*mm)}
};

// ----------------------------------------------------------------------------
//  6. Alpha source (GDML) placement inside the LXe
// ----------------------------------------------------------------------------
inline constexpr G4double sourceX = GEO_V(0.*mm);
inline constexpr G4double sourceY = GEO_V(0.*mm);
inline constexpr G4double sourceZ = GEO_V(-80.*mm);

// ----------------------------------------------------------------------------
//  7. Environment and numerics
// ----------------------------------------------------------------------------
// Material on the beamline side of the Al window.  CONFIRM: vacuum or air?
inline constexpr const char* worldMaterial = "G4_Galactic";
inline constexpr G4double worldMargin = 200.*mm;
// Production cut inside the foils (must be << foil thickness; the energy
// threshold is floored by Geant4 at ~1 keV regardless).
inline constexpr G4double windowProductionCut = 2.*um;

#undef GEO_V

// ============================================================================
//  Derived quantities and helpers (no CAD numbers below this line)
// ============================================================================
inline G4double outerR() { return outer.cavityR + outer.wallT; }  // outer surface of outer wall
inline G4double innerR() { return inner.cavityR + inner.wallT; }  // outer surface of inner wall

inline G4double worldHalfX() { return std::max(outerR(), outerPort.flangeS1) + worldMargin; }
inline G4double worldHalfY() { return std::max(outerR(), outerPort.flangeOR) + 100.*mm; }
inline G4double worldHalfZ()
{
  return std::max(std::abs(outer.cavityZmin - outer.bottomT),
                  std::abs(outer.cavityZmax + outer.topT)) + 100.*mm;
}
// Upstream-most material is the outer flange face; start beams/scans before it.
inline G4double suggestedBeamStartX() { return -(outerPort.flangeS1 + 50.*mm); }

inline G4RotationMatrix PlaneRotation(const PlaneSpec& p)
{
  G4RotationMatrix r;
  r.rotateX(p.rotXdeg * deg);
  r.rotateY(p.rotYdeg * deg);
  r.rotateZ(p.rotZdeg * deg);
  return r;
}

// Returns human-readable problems; empty == geometry is complete and consistent.
inline std::vector<std::string> Validate()
{
  std::vector<std::string> bad;

  auto fin = [&](const std::string& n, G4double v) {
    if (!std::isfinite(v)) bad.push_back(n + " is not set (still TBD)");
  };
  auto finVessel = [&](const std::string& n, const Vessel& v) {
    fin(n + ".cavityR", v.cavityR);       fin(n + ".cavityZmin", v.cavityZmin);
    fin(n + ".cavityZmax", v.cavityZmax); fin(n + ".wallT", v.wallT);
    fin(n + ".bottomT", v.bottomT);       fin(n + ".topT", v.topT);
  };
  auto finPort = [&](const std::string& n, const Port& p) {
    fin(n + ".boreR", p.boreR);     fin(n + ".stubOR", p.stubOR);
    fin(n + ".flangeOR", p.flangeOR);
    fin(n + ".flangeS0", p.flangeS0); fin(n + ".flangeS1", p.flangeS1);
  };
  finVessel("outer", outer);
  finVessel("inner", inner);
  finPort("outerPort", outerPort);
  finPort("innerPort", innerPort);
  fin("alWindow.thickness", alWindow.thickness); fin("alWindow.outerFaceS", alWindow.outerFaceS);
  fin("tiWindow.thickness", tiWindow.thickness); fin("tiWindow.outerFaceS", tiWindow.outerFaceS);
  if (sleeve.use) {
    fin("sleeve.IR", sleeve.IR); fin("sleeve.OR", sleeve.OR);
    fin("sleeve.s0", sleeve.s0); fin("sleeve.s1", sleeve.s1);
  }
  for (G4int i = 0; i < nPlanes; ++i) {
    const auto& p = planes[i];
    const std::string n = std::string("planes[") + std::to_string(i) + "] (" + p.name + ")";
    fin(n + ".cx", p.cx); fin(n + ".cy", p.cy); fin(n + ".cz", p.cz);
    fin(n + ".rotX", p.rotXdeg); fin(n + ".rotY", p.rotYdeg); fin(n + ".rotZ", p.rotZdeg);
    fin(n + ".halfU", p.halfU); fin(n + ".halfV", p.halfV); fin(n + ".thickness", p.thickness);
  }
  fin("sourceX", sourceX); fin("sourceY", sourceY); fin("sourceZ", sourceZ);
  if (!bad.empty()) return bad;           // geometric checks are meaningless with NaNs

  auto req = [&](G4bool ok, const std::string& msg) { if (!ok) bad.push_back(msg); };

  // --- basic sanity -------------------------------------------------------
  for (const Vessel* v : {&outer, &inner}) {
    const std::string n = (v == &outer) ? "outer" : "inner";
    req(v->cavityR > 0 && v->wallT > 0 && v->bottomT > 0 && v->topT > 0,
        n + ": radii/thicknesses must be > 0");
    req(v->cavityZmax > v->cavityZmin, n + ": cavityZmax must exceed cavityZmin");
  }
  for (const Port* p : {&outerPort, &innerPort}) {
    const std::string n = (p == &outerPort) ? "outerPort" : "innerPort";
    req(p->boreR > 0 && p->boreR < p->stubOR && p->stubOR < p->flangeOR,
        n + ": need 0 < boreR < stubOR < flangeOR");
    req(p->flangeS1 > p->flangeS0, n + ": flangeS1 must exceed flangeS0");
  }

  // --- nesting ------------------------------------------------------------
  req(innerR() + 1.*mm <= outer.cavityR, "inner vessel (incl. wall) must fit inside the outer cavity with >= 1 mm clearance");
  req(inner.cavityZmin - inner.bottomT > outer.cavityZmin, "inner vessel bottom plate collides with the outer cavity");
  req(inner.cavityZmax + inner.topT    < outer.cavityZmax, "inner vessel top plate collides with the outer cavity");

  // --- port geometry ------------------------------------------------------
  // Flange inner face must clear the vessel's outer cylinder over its whole disc.
  req(outerPort.flangeS0 >= outerR(), "outerPort.flangeS0 must be >= outer vessel outer radius");
  req(innerPort.flangeS0 >= innerR(), "innerPort.flangeS0 must be >= inner vessel outer radius");
  req(outerPort.stubOR < outerR() && innerPort.stubOR < innerR(), "stub OD larger than the vessel it is welded to");
  // Inner flange must sit inside the outer cavity (curved wall, hence the sqrt).
  if (innerPort.flangeOR < outer.cavityR)
    req(innerPort.flangeS1 + 1.*mm <= std::sqrt(outer.cavityR*outer.cavityR - innerPort.flangeOR*innerPort.flangeOR),
        "inner flange (outer face + disc radius) pokes into the outer vessel wall");
  else
    bad.push_back("innerPort.flangeOR >= outer.cavityR");
  // Port is at z = 0: its disc must lie within the straight part of each vessel.
  req(-innerPort.stubOR > inner.cavityZmin && innerPort.stubOR < inner.cavityZmax, "inner stub spans past the inner vessel's z extent");
  req(-outerPort.flangeOR > outer.cavityZmin && outerPort.flangeOR < outer.cavityZmax, "outer flange spans past the outer vessel's z extent");
  req(-innerPort.flangeOR > outer.cavityZmin && innerPort.flangeOR < outer.cavityZmax, "inner flange spans past the outer vessel's z extent");

  // --- windows live inside their flange bore ------------------------------
  req(alWindow.thickness > 0 && alWindow.outerFaceS <= outerPort.flangeS1 &&
      alWindow.outerFaceS - alWindow.thickness >= outerPort.flangeS0,
      "alWindow must lie within [outerPort.flangeS0, outerPort.flangeS1]");
  req(tiWindow.thickness > 0 && tiWindow.outerFaceS <= innerPort.flangeS1 &&
      tiWindow.outerFaceS - tiWindow.thickness >= innerPort.flangeS0,
      "tiWindow must lie within [innerPort.flangeS0, innerPort.flangeS1]");

  // --- sleeve -------------------------------------------------------------
  if (sleeve.use) {
    req(sleeve.IR > 0 && sleeve.IR < sleeve.OR, "sleeve: need 0 < IR < OR");
    req(sleeve.OR < outerPort.boreR, "sleeve.OR must be < outerPort.boreR (it sits inside the outer stub)");
    req(sleeve.s1 > sleeve.s0, "sleeve: s1 must exceed s0");
    req(sleeve.s0 >= innerPort.flangeS1, "sleeve overlaps the inner flange (s0 < innerPort.flangeS1)");
    req(sleeve.s1 <= alWindow.outerFaceS - alWindow.thickness, "sleeve runs into the Al window");
  }

  // --- planes ---------------------------------------------------------------
  for (G4int i = 0; i < nPlanes; ++i) {
    const auto& p = planes[i];
    req(p.halfU > 0 && p.halfV > 0 && p.thickness > 0,
        std::string("planes[") + std::to_string(i) + "]: halfU, halfV, thickness must be > 0");
  }
  return bad;
}

// Written next to the output at construction so every ROOT file can be tied to
// the exact geometry that produced it (analysis scripts read this).
inline std::string ManifestJSON()
{
  std::ostringstream o;
  o.precision(10);
  auto vessel = [&](const char* n, const Vessel& v) {
    o << "  \"" << n << "\": {\"cavityR\": " << v.cavityR << ", \"cavityZmin\": " << v.cavityZmin
      << ", \"cavityZmax\": " << v.cavityZmax << ", \"wallT\": " << v.wallT << ", \"bottomT\": " << v.bottomT
      << ", \"topT\": " << v.topT << ", \"material\": \"" << v.material << "\"},\n";
  };
  auto port = [&](const char* n, const Port& p) {
    o << "  \"" << n << "\": {\"boreR\": " << p.boreR << ", \"stubOR\": " << p.stubOR << ", \"flangeOR\": " << p.flangeOR
      << ", \"flangeS0\": " << p.flangeS0 << ", \"flangeS1\": " << p.flangeS1 << ", \"material\": \"" << p.material << "\"},\n";
  };
  auto window = [&](const char* n, const Window& w) {
    o << "  \"" << n << "\": {\"thickness\": " << w.thickness << ", \"outerFaceS\": " << w.outerFaceS
      << ", \"material\": \"" << w.material << "\"},\n";
  };
  o << "{\n  \"isDemo\": " << (isDemo ? "true" : "false") << ",\n  \"units\": \"mm\",\n";
  vessel("outer", outer);
  vessel("inner", inner);
  port("outerPort", outerPort);
  port("innerPort", innerPort);
  window("alWindow", alWindow);
  window("tiWindow", tiWindow);
  o << "  \"sleeve\": {\"use\": " << (sleeve.use ? "true" : "false") << ", \"IR\": " << sleeve.IR << ", \"OR\": " << sleeve.OR
    << ", \"s0\": " << sleeve.s0 << ", \"s1\": " << sleeve.s1 << ", \"material\": \"" << sleeve.material << "\"},\n";
  o << "  \"worldMaterial\": \"" << worldMaterial << "\",\n";
  o << "  \"suggestedBeamStartX\": " << suggestedBeamStartX() << ",\n";
  o << "  \"planes\": [\n";
  for (G4int i = 0; i < nPlanes; ++i) {
    const auto& p = planes[i];
    const G4ThreeVector n = PlaneRotation(p) * G4ThreeVector(0, 0, 1);
    o << "    {\"id\": " << i << ", \"name\": \"" << p.name << "\", \"centre\": [" << p.cx << ", " << p.cy << ", " << p.cz
      << "], \"normal\": [" << n.x() << ", " << n.y() << ", " << n.z() << "], \"halfU\": " << p.halfU
      << ", \"halfV\": " << p.halfV << ", \"thickness\": " << p.thickness << "}" << (i + 1 < nPlanes ? "," : "") << "\n";
  }
  o << "  ]\n}\n";
  return o.str();
}

}  // namespace Geo

#endif
