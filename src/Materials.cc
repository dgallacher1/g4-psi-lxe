#include "Materials.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Exception.hh"
#include "G4SystemOfUnits.hh"

#include <string>
#include <utility>
#include <vector>

namespace Materials {
namespace {

// A material is a density and a list of (NIST element symbol, mass fraction). Fractions are normalised
// when the material is built, so a list that does not add up to exactly 1 is still consistent; the
// "balance" element (listed first, fraction ignored) takes whatever the others leave.
struct Component { const char* element; double massFraction; };

struct Recipe {
  const char* name;
  double density;                       // internal units
  bool placeholder;                     // composition or density not confirmed yet
  const char* note;
  std::vector<Component> components;    // first entry = balance (its fraction is computed)
};

// ---------------------------------------------------------------------------------------------
//  EDIT HERE
// ---------------------------------------------------------------------------------------------
const std::vector<Recipe>& Recipes()
{
  static const std::vector<Recipe> r = {

    // AISI / ASTM 304 stainless steel, nominal composition (mass %): Cr 18-20, Ni 8-10.5, Mn <= 2, Si <= 0.75,
    // C <= 0.08, P <= 0.045, S <= 0.03, balance Fe. The values below are mid-range typicals; the density is the
    // usual 7.93 g/cm3 (G4_STAINLESS-STEEL is Fe 74 / Cr 18 / Ni 8 at 8.0 g/cm3).
    {"SS304", 7.93 * g / cm3, false, "nominal 304",
     {{"Fe", 0.0}, {"Cr", 0.1900}, {"Ni", 0.0925}, {"Mn", 0.0100}, {"Si", 0.0050},
      {"C", 0.0004}, {"P", 0.0002}, {"S", 0.0001}}},

    // Ti window alloy: Ti-6Al-4V (mass fractions Ti 0.90, Al 0.06, V 0.04, as in the GPU light simulation).
    // Density 4.5 g/cm3, the value used in the GPU light simulation.
    {"TiWindowAlloy", 4.5 * g / cm3, false, "Ti-6Al-4V, density 4.5 g/cm3 (as in the GPU simulation)",
     {{"Ti", 0.0}, {"Al", 0.06}, {"V", 0.04}}},
  };
  return r;
}
// ---------------------------------------------------------------------------------------------

void Build(const Recipe& r)
{
  if (G4Material::GetMaterial(r.name, false)) return;
  auto* nist = G4NistManager::Instance();

  double others = 0.;
  for (std::size_t i = 1; i < r.components.size(); ++i) others += r.components[i].massFraction;
  if (r.components.empty() || others >= 1.)
    G4Exception("Materials::Define", "Mat001", FatalException,
                (std::string("material ") + r.name + ": the non-balance fractions must add up to less than 1").c_str());

  auto* m = new G4Material(r.name, r.density, static_cast<G4int>(r.components.size()), kStateSolid);
  for (std::size_t i = 0; i < r.components.size(); ++i) {
    auto* el = nist->FindOrBuildElement(r.components[i].element);
    if (!el)
      G4Exception("Materials::Define", "Mat002", FatalException,
                  (std::string("material ") + r.name + ": unknown element " + r.components[i].element).c_str());
    m->AddElement(el, i == 0 ? 1. - others : r.components[i].massFraction);
  }

  G4cout << "Material " << r.name << ": " << r.density / (g / cm3) << " g/cm3, " << r.components.size()
         << " components" << (r.placeholder ? "   *** PLACEHOLDER: " : "   (") << r.note
         << (r.placeholder ? " ***" : ")") << G4endl;
}

}  // namespace

void Define()
{
  for (const auto& r : Recipes()) Build(r);
}

bool IsPlaceholder(const std::string& materialName)
{
  for (const auto& r : Recipes())
    if (materialName == r.name) return r.placeholder;
  return false;
}

}  // namespace Materials
