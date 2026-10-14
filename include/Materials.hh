#ifndef Materials_hh
#define Materials_hh 1

#include "globals.hh"
#include <string>

// Custom materials that are not in the NIST database with the right composition.
// Defined once on the master before the geometry is built (xenonSim.cc); the geometry then asks for them
// by name, like any NIST material: GeometryParameters.hh uses "SS304" and "TiWindowAlloy".
//
// To change a composition or a density edit the tables in Materials.cc, nothing else.
namespace Materials {

// Idempotent: creates SS304 and TiWindowAlloy if they do not exist yet.
void Define();

// true while the named custom material is a placeholder (composition not confirmed): printed at start-up
// and recorded in the run_info ntuple, so a placeholder can never be mistaken for the final alloy.
bool IsPlaceholder(const std::string& materialName);

}  // namespace Materials

#endif
