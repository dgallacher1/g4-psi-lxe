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
//        outer CF flange / Al window  ->  insulating vacuum  ->  inner CF flange /
//        Ti window  ->  LXe.
//    CAD positions on the port are measured from the vessel wall; convert to s
//    with s = (vessel outer radius) + (CAD value from the wall OD), or
//    s = (vessel inner radius) + (CAD value from the wall ID).
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
#include "G4PhysicalConstants.hh"
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
  G4double halfU, halfV;              // circular plane: halfU = halfV = radius
  G4double thickness;                 // silicon thickness (a plane is a real, interacting slab)
  G4bool   circular;                  // true: disc of radius halfU (G4Tubs); false: rectangle (G4Box)
};

inline G4double PlaneArea(const PlaneSpec& p)
{
  return p.circular ? CLHEP::pi * p.halfU * p.halfU : 4. * p.halfU * p.halfV;
}

// ----------------------------------------------------------------------------
//  1. Vessels   (CAD: cut-away section, vertical cylinders)
// ----------------------------------------------------------------------------
//  CAD values (Fusion).  Radii and z extents as measured.
inline constexpr Vessel outer = {
  /*cavityR   */ 82.5*mm,
  /*cavityZmin*/ -187.*mm,
  /*cavityZmax*/ 187.*mm,
  /*wallT     */ 6.*mm,
  /*bottomT   */ 4.*mm,
  /*topT      */ 4.*mm,
  "G4_STAINLESS-STEEL"
};

inline constexpr Vessel inner = {
  /*cavityR   */ 50.*mm,
  /*cavityZmin*/ -141.2*mm,
  /*cavityZmax*/ 141.2*mm,
  /*wallT     */ 2.*mm,
  /*bottomT   */ 4.*mm,
  /*topT      */ 4.*mm,
  "G4_STAINLESS-STEEL"
};

// ----------------------------------------------------------------------------
//  2. Beam port on the OUTER vessel  (DN40CF-style flange)
//     stub is welded to the wall; the flange is a simple additive ring.
//     CAD s values are measured from the outer vessel OD (88.5 mm):
//       flange 50.5 / 74.6 -> s = 139.0 / 163.1,  Al window 62.63 -> s = 151.13
//     The "flange" spans two half-flanges; the Al foil is clamped between them,
//     so it sits in the middle of the ring.  Knife edge, gasket groove and bolt
//     holes are not modelled (they are away from the beam axis).
// ----------------------------------------------------------------------------
inline constexpr Port outerPort = {
  /*boreR   */ 19.2*mm,
  /*stubOR  */ 21.2*mm,                 // CAD gives the tube OD = 42.4 mm
  /*flangeOR*/ 35.0*mm,
  /*flangeS0*/ 139.0*mm,
  /*flangeS1*/ 163.1*mm,
  "G4_STAINLESS-STEEL"
};

// First (thin) window: aluminium, separates the jacket vacuum from the beamline side.
inline constexpr Window alWindow = {
  /*thickness */ 0.1*mm,
  /*outerFaceS*/ 151.13*mm,
  "G4_Al"
};

// ----------------------------------------------------------------------------
//  3. Sleeve: not present in the CAD (the dark tube in the section view is the
//     inside surface of the outer stub).  Kept as an option; leave use = false.
// ----------------------------------------------------------------------------
inline constexpr Sleeve sleeve = {
  /*use*/ false,
  /*IR */ 0.*mm,
  /*OR */ 0.*mm,
  /*s0 */ 0.*mm,
  /*s1 */ 0.*mm,
  "G4_Al"
};

// ----------------------------------------------------------------------------
//  4. Port on the INNER vessel: non-standard CF, concentric with the outer
//     port, NOT attached to it (cryo insulation gap).
//
//     There is no separate stub.  The flange fills an opening in the wall (opening
//     radius = flangeOR) and reaches 3.77 mm into the cavity.  From the cavity
//     side it has a POCKET of radius innerPocketR that holds the SiPM board; the
//     pocket floor is the Ti foil's inner face, and the foil is welded across the
//     narrower boreR beyond it.  The LXe fills the wall opening and the pocket
//     (it is built minus the flange metal), so the SiPM and the foil see LXe.
//
//       s (from the CAD)                                        value used
//       flange inner face  -3.77 from wall ID (50)       ->      46.23
//         (cylindrical, concentric with the vessel: s on the port axis = its radius)
//       flange outer face  16.71 from wall OD (52)       ->      68.71
//       Ti outer face       8.50 from wall OD (52)       ->      60.50
//
//     stubOR is unused for this port (the CAD "stub OR 20 mm" is the pocket).
// ----------------------------------------------------------------------------
inline constexpr Port innerPort = {
  /*boreR   */ 12.5*mm,
  /*stubOR  */ 0.*mm,                   // unused (see above)
  /*flangeOR*/ 34.8*mm,
  /*flangeS0*/ 46.23*mm,
  /*flangeS1*/ 68.71*mm,
  "G4_STAINLESS-STEEL"
};
inline constexpr G4double innerPocketR = 20.*mm;   // cavity-side opening of the inner flange

inline constexpr Window tiWindow = {
  /*thickness */ 0.5*mm,               // 500 um
  /*outerFaceS*/ 60.5*mm,
  "G4_Ti"
};

// ----------------------------------------------------------------------------
//  5. Fluence planes standing in for the SiPMs (silicon slabs, children of the LXe)
//     They are real silicon, so they scatter, absorb and shadow like a SiPM would.
// ----------------------------------------------------------------------------
inline constexpr const char* planeMaterial = "G4_Si";
inline constexpr G4int nPlanes = 2;

//  CAD positions use Fusion's X, which points TOWARD the port; ours points away
//  from it, so x_here = -X_cad (a 180 deg turn about z; all CAD y are 0).
inline constexpr PlaneSpec planes[nPlanes] = {
  // name    cx             cy           cz             rotX         rotY         rotZ        halfU          halfV          thick (mm)       circular
  // SiPM1: disc r = 18 mm in the flange pocket, 4 mm behind the Ti foil, normal toward the port (-x).
  //        Thickness is NOT from the CAD: 0.05 mm silicon assumed, replace with the real die thickness.
  {"SiPM1", -56.41*mm,     0.*mm,        0.*mm,         0.,          -90.,        0.,         18.*mm,        18.*mm,        0.05*mm,         true},
  // SiPM2: rectangle 55 x 73 mm, normal +z.  Its far edge (the side away from the port) is 13.80 mm from the
  //        inner wall ID (x = +50), so the far edge is at x = 36.2 and the centre at 36.2 - 27.5 = 8.7.
  //        The real board is not rectangular: the length (y) is cut from 73 to 63 mm (halfV 31.5) so the
  //        corners fit in the cylinder (r = 48.0 mm); the SiPMs are inside this area.
  //        z = -32.63 is the CAD board plane; the SiPM silicon may sit higher, check.  Thickness assumed.
  {"SiPM2", 8.7*mm,        0.*mm,       -32.63*mm,      0.,          0.,          0.,         27.5*mm,       31.5*mm,       0.05*mm,         false}
};

// ----------------------------------------------------------------------------
//  6. Alpha source (GDML) placement inside the LXe   (CAD X = +30.31 -> x = -30.31)
// ----------------------------------------------------------------------------
inline constexpr G4double sourceX = -30.31*mm;
inline constexpr G4double sourceY = 0.*mm;
inline constexpr G4double sourceZ = 0.*mm;

// ----------------------------------------------------------------------------
//  7. Environment and numerics
// ----------------------------------------------------------------------------
// Material outside the Al window (and around the vessels): air (confirmed).
// The air path between the beam-profile monitor and the window scatters the
// beam, so the beam start plane (macros/beam.mac, /beam/offsetX) should sit at
// the monitor's position, not at an arbitrary distance upstream.
inline constexpr const char* worldMaterial = "G4_AIR";
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
  req(outerPort.boreR > 0 && outerPort.boreR < outerPort.stubOR && outerPort.stubOR < outerPort.flangeOR,
      "outerPort: need 0 < boreR < stubOR < flangeOR");
  req(innerPort.boreR > 0 && innerPort.boreR < innerPocketR && innerPocketR < innerPort.flangeOR,
      "innerPort: need 0 < boreR < innerPocketR < flangeOR");
  for (const Port* p : {&outerPort, &innerPort})
    req(p->flangeS1 > p->flangeS0, (p == &outerPort ? "outerPort" : "innerPort") + std::string(": flangeS1 must exceed flangeS0"));

  // --- nesting ------------------------------------------------------------
  req(innerR() + 1.*mm <= outer.cavityR, "inner vessel (incl. wall) must fit inside the outer cavity with >= 1 mm clearance");
  req(inner.cavityZmin - inner.bottomT > outer.cavityZmin, "inner vessel bottom plate collides with the outer cavity");
  req(inner.cavityZmax + inner.topT    < outer.cavityZmax, "inner vessel top plate collides with the outer cavity");

  // --- port geometry ------------------------------------------------------
  // Flange inner face must clear the vessel's outer cylinder over its whole disc.
  req(outerPort.flangeS0 >= outerR(), "outerPort.flangeS0 must be >= outer vessel outer radius");
  // (The inner flange may reach into the cavity: the LXe is built minus the flange.)
  req(innerPort.flangeS0 > innerPort.flangeOR, "innerPort.flangeS0 (radius of the curved inner face) must exceed flangeOR");
  req(outerPort.stubOR < outerR(), "outer stub OD larger than the vessel it is welded to");
  // Inner flange must sit inside the outer cavity (curved wall, hence the sqrt).
  if (innerPort.flangeOR < outer.cavityR)
    req(innerPort.flangeS1 + 1.*mm <= std::sqrt(outer.cavityR*outer.cavityR - innerPort.flangeOR*innerPort.flangeOR),
        "inner flange (outer face + disc radius) pokes into the outer vessel wall");
  else
    bad.push_back("innerPort.flangeOR >= outer.cavityR");
  // Port is at z = 0: its disc must lie within the straight part of each vessel.
  req(-innerPort.flangeOR > inner.cavityZmin && innerPort.flangeOR < inner.cavityZmax, "inner flange spans past the inner vessel's z extent");
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
    req(!p.circular || p.halfU == p.halfV,
        std::string("planes[") + std::to_string(i) + "]: a circular plane needs halfU == halfV (= radius)");
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
  o << "  \"innerPocketR\": " << innerPocketR << ",\n";
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
      << ", \"halfV\": " << p.halfV << ", \"thickness\": " << p.thickness
      << ", \"circular\": " << (p.circular ? "true" : "false") << ", \"area\": " << PlaneArea(p)
      << "}" << (i + 1 < nPlanes ? "," : "") << "\n";
  }
  o << "  ]\n}\n";
  return o.str();
}

}  // namespace Geo

#endif
