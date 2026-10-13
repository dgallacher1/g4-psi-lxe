// geometry_probe.cc -- checks the solids built from GeometryParameters.hh WITHOUT
// running a simulation.  Uses exactly the same builders as the application
// (GeometryBuilders.hh), with probe points derived from the parameters, so it
// stays valid after you replace the placeholders with CAD numbers.
//
//   ./geometryProbe          -> exit code 0 if everything passes
//
// What it checks
//   * Validate(): header is complete and self-consistent
//   * orientation: port parts are on the -x side only (beam enters from -x)
//   * booleans: port bore cut through each wall, stub root cut by the vessel
//   * foils: on-axis thickness along +x equals the configured thickness
//   * LXe / jacket vacuum extend through the bores up to their foils
//   * fluence planes: normals as documented, and fully inside the LXe
#include "GeometryParameters.hh"
#include "GeometryBuilders.hh"

#include "G4VSolid.hh"
#include <cmath>
#include <iostream>
#include <string>

namespace {

int nFail = 0, nPass = 0;

void check(bool ok, const char* what)
{
  (ok ? nPass : nFail)++;
  std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << "\n";
}

// Point on the port at distance s from the vessel axis, offset (dy, dz) from the port axis.
G4ThreeVector P(G4double s, G4double dy = 0., G4double dz = 0.) { return {-s, dy, dz}; }
G4ThreeVector Pmirror(G4double s, G4double dy = 0., G4double dz = 0.) { return {+s, dy, dz}; }

bool in(const G4VSolid* v, const G4ThreeVector& p)  { return v->Inside(p) == kInside; }
bool out(const G4VSolid* v, const G4ThreeVector& p) { return v->Inside(p) == kOutside; }

// Fire a ray along +x on the port axis from 1 mm upstream of the foil's outer face.
//   entryOffset : (distance to the first surface) - 1 mm  -> 0 if the face is at the configured s
//   thickness   : path length inside the solid
struct Ray { G4double entryOffset, thickness; };
Ray AxisRay(G4VSolid* v, G4double outerFaceS)
{
  const G4ThreeVector dir(1., 0., 0.), p0 = P(outerFaceS + 1.*mm);
  const G4double eps = 1e-6*mm;
  const G4double d = v->DistanceToIn(p0, dir);
  const G4ThreeVector p1 = p0 + (d + eps) * dir;        // just inside the solid
  return {d - 1.*mm, v->DistanceToOut(p1, dir) + eps};
}

}  // namespace

int main(int argc, char** argv)
{
  using namespace Geo;

  if (argc > 1 && std::string(argv[1]) == "--manifest") {   // dump geometry_manifest.json to stdout
    if (!Validate().empty()) { std::cerr << "header incomplete\n"; return 1; }
    std::cout << ManifestJSON();
    return 0;
  }

  std::cout << "geometry_probe (" << (isDemo ? "DEMO placeholders" : "CAD values") << ")\n\n[Validate]\n";
  const auto problems = Validate();
  for (const auto& p : problems) std::cout << "  problem: " << p << "\n";
  check(problems.empty(), "header complete and self-consistent");
  if (!problems.empty()) { std::cout << "\nstopping: fix the header first\n"; return 1; }

  const auto S = GeoBuild::Build();
  const G4double RO = outerR(), RI = innerR();
  const auto& OP = outerPort; const auto& IP = innerPort;
  const G4double yO = 0.5 * (OP.boreR + OP.stubOR);          // mid-annulus of the outer stub
  const G4double fO = 0.5 * (OP.boreR + OP.flangeOR);        // mid-annulus of the flange rings
  const G4double fI = 0.5 * (innerPocketR + IP.flangeOR);    // mid-annulus of the inner flange's pocket ring
  const G4double sFO = 0.5 * (OP.flangeS0 + OP.flangeS1);    // flange mid-plane
  const G4double sFI = 0.5 * (IP.flangeS0 + IP.flangeS1);

  std::cout << "\n[Orientation: port is on -x only]\n";
  check(in(S.outerFlange, P(sFO, fO)) && out(S.outerFlange, Pmirror(sFO, fO)), "outer flange at -x, not at +x");
  check(in(S.innerFlange, P(sFI, fI)) && out(S.innerFlange, Pmirror(sFI, fI)), "inner flange at -x, not at +x");
  check(in(S.outerFlange, P(sFO, 0., fO)), "outer flange symmetric in z (ring, not slab)");
  check(out(S.outerFlange, P(sFO, 0.5 * OP.boreR)), "outer flange has its bore");

  std::cout << "\n[Wall bores]\n";
  {
    const G4double rw = outer.cavityR + 0.5 * outer.wallT;
    const G4double yIn = OP.boreR - 0.1*mm, yOut = OP.boreR + 0.1*mm;
    const auto wallPt = [&](G4double y) { return G4ThreeVector(-std::sqrt(rw*rw - y*y), y, 0.); };
    check(out(S.outerWall, wallPt(yIn)),  "outer wall open inside the port bore (-x)");
    check(in(S.outerWall,  wallPt(yOut)), "outer wall intact just outside the bore radius");
    check(in(S.outerWall,  G4ThreeVector(+rw, 0., 0.)), "outer wall intact at +x (no bore on the far side)");
  }
  {
    // The inner wall's opening is as wide as the flange that fills it.
    const G4double rw = inner.cavityR + 0.5 * inner.wallT;
    const G4double yIn = IP.flangeOR - 0.1*mm, yOut = IP.flangeOR + 0.1*mm;
    const auto wallPt = [&](G4double y) { return G4ThreeVector(-std::sqrt(rw*rw - y*y), y, 0.); };
    check(out(S.innerWall, wallPt(yIn)),  "inner wall open inside the flange opening (-x)");
    check(in(S.innerWall,  wallPt(yOut)), "inner wall intact just outside the flange opening");
    check(in(S.innerWall,  G4ThreeVector(+rw, 0., 0.)), "inner wall intact at +x");
  }

  std::cout << "\n[Stub roots cut by the vessel outer cylinder]\n";
  {
    const G4double sEdge = std::sqrt(RO*RO - yO*yO);
    check(out(S.outerStub, P(sEdge - 0.5*mm, yO)), "outer stub removed inside the vessel outer cylinder");
    check(in(S.outerStub,  P(sEdge + 2.*mm, yO)), "outer stub present just outside the vessel");
    check(in(S.outerStub,  P(OP.flangeS0 - 1.*mm, yO)), "outer stub reaches the flange");
    check(out(S.outerStub, Pmirror(sEdge + 2.*mm, yO)), "outer stub not mirrored to +x");
    check(out(S.outerWall, P(sEdge + 2.*mm, yO)), "stub region is not also wall (no double material)");
  }
  check(S.innerStub == nullptr, "inner port has no stub (flange fills the wall opening)");

  std::cout << "\n[Foils]\n";
  {
    const Ray rAl = AxisRay(S.alWindow, alWindow.outerFaceS);
    const Ray rTi = AxisRay(S.tiWindow, tiWindow.outerFaceS);
    std::cout << "        on-axis thickness  Al = " << rAl.thickness / um << " um (configured " << alWindow.thickness / um
              << ")   Ti = " << rTi.thickness / um << " um (configured " << tiWindow.thickness / um << ")\n";
    check(std::abs(rAl.entryOffset) < 1e-3*mm, "Al foil outer face is at alWindow.outerFaceS");
    check(std::abs(rTi.entryOffset) < 1e-3*mm, "Ti foil outer face is at tiWindow.outerFaceS");
    check(std::abs(rAl.thickness - alWindow.thickness) < 1e-3*mm, "Al foil on-axis thickness matches the header");
    check(std::abs(rTi.thickness - tiWindow.thickness) < 1e-3*mm, "Ti foil on-axis thickness matches the header");
    check(out(S.alWindow, P(alWindow.outerFaceS - 0.5*alWindow.thickness, OP.boreR + 0.1*mm)), "Al foil limited to the bore radius");
    check(out(S.tiWindow, P(tiWindow.outerFaceS - 0.5*tiWindow.thickness, IP.boreR + 0.1*mm)), "Ti foil limited to the bore radius");
    check(out(S.alWindow, Pmirror(alWindow.outerFaceS - 0.5*alWindow.thickness)), "no Al foil at +x");
  }

  std::cout << "\n[Vacuum and LXe reach their foils]\n";
  check(in(S.jacket, G4ThreeVector(0., 0., 0.)),                       "jacket vacuum fills the outer cavity");
  check(in(S.jacket, P(0.5 * (OP.flangeS0 + S.sAlIn))),               "jacket vacuum reaches through the outer port bore");
  check(out(S.jacket, P(alWindow.outerFaceS + 0.5*mm)),               "jacket vacuum stops at the Al foil");
  check(out(S.jacket, P(0.5 * (RO + OP.flangeS0), OP.boreR + 0.1*mm)), "jacket vacuum confined to the bore radius");
  check(in(S.lxe, G4ThreeVector(0., 0., 0.)),                          "LXe fills the inner cavity");
  check(in(S.lxe, P(0.5 * (inner.cavityR + S.sTiIn))),                "LXe reaches through the inner port bore");
  check(out(S.lxe, P(S.sTiIn + 0.5 * tiWindow.thickness)),            "LXe stops at the Ti foil");
  check(in(S.lxe,  P(S.sTiIn - 1.*mm, innerPocketR - 0.1*mm)),        "LXe fills the flange pocket");
  check(out(S.lxe, P(S.sTiIn - 1.*mm, innerPocketR + 0.1*mm)),        "LXe stops at the pocket wall (flange metal)");
  check(out(S.lxe, P(S.sTiIn + 1.*mm, IP.boreR - 0.1*mm)),            "no LXe beyond the Ti foil");
  {
    // the flange's cavity-side face is a cylinder of radius flangeS0 about the vessel axis
    const G4double sFace = std::sqrt(IP.flangeS0*IP.flangeS0 - fI*fI);       // face s at offset y = fI
    const G4double sEdge = std::sqrt(IP.flangeS0*IP.flangeS0 - (IP.flangeOR - 0.5*mm)*(IP.flangeOR - 0.5*mm));
    check(in(S.lxe,  P(sFace - 0.5*mm, fI)) && out(S.lxe, P(sFace + 0.5*mm, fI)),
          "inner flange face is curved with the vessel (LXe stops at it, mid-annulus)");
    check(in(S.lxe,  P(sEdge - 0.5*mm, IP.flangeOR - 0.5*mm)) && out(S.lxe, P(sEdge + 0.5*mm, IP.flangeOR - 0.5*mm)),
          "inner flange reaches the cavity at the edge of the wall opening (no gap)");
    check(in(S.innerFlange, P(sEdge + 0.5*mm, IP.flangeOR - 0.5*mm)), "inner flange fills the wall opening out to its edge");
  }
  check(out(S.lxe, G4ThreeVector(+inner.cavityR + 0.5*mm, 0., 0.)),   "LXe does not leak through the +x wall");
  if (S.sleeve) {
    const G4double sm = 0.5 * (sleeve.s0 + sleeve.s1), rm = 0.5 * (sleeve.IR + sleeve.OR);
    check(in(S.sleeve, P(sm, rm)) && out(S.sleeve, P(sm, 0.5 * sleeve.IR)), "sleeve is a hollow tube on -x");
    check(out(S.sleeve, Pmirror(sm, rm)), "no sleeve at +x");
  }

  std::cout << "\n[Fluence planes]\n";
  for (G4int i = 0; i < nPlanes; ++i) {
    const auto& p = planes[i];
    const G4RotationMatrix R = PlaneRotation(p);
    const G4ThreeVector n = R * G4ThreeVector(0, 0, 1), U = R * G4ThreeVector(1, 0, 0), V = R * G4ThreeVector(0, 1, 0);
    std::cout << "        " << p.name << ": normal (" << n.x() << ", " << n.y() << ", " << n.z() << ")  U ("
              << U.x() << ", " << U.y() << ", " << U.z() << ")  V (" << V.x() << ", " << V.y() << ", " << V.z() << ")\n";
    bool allIn = true;
    // Rectangle: the 8 corners.  Disc: the 4 rim points on U and V, on both faces.
    for (int a : {-1, 0, 1}) for (int b : {-1, 0, 1}) for (int c : {-1, 1}) {
      if (p.circular ? (a * b != 0 || (a == 0 && b == 0)) : (a == 0 || b == 0)) continue;
      const G4ThreeVector corner = R * G4ThreeVector(a * p.halfU, b * p.halfV, c * 0.5 * p.thickness)
                                 + G4ThreeVector(p.cx, p.cy, p.cz);
      if (S.lxe->Inside(corner) == kOutside) allIn = false;
    }
    check(allIn, (std::string(p.name) + " lies entirely inside the LXe").c_str());
  }

  std::cout << "\n" << nPass << " passed, " << nFail << " failed\n";
  return nFail == 0 ? 0 : 1;
}
