#ifndef OutputSchema_hh
#define OutputSchema_hh 1

// ---------------------------------------------------------------------------------------------
//  THE output format: every ntuple and every column is defined here, once.
//
//  * RunAction creates the ntuples from these tables (Out::Create<>).
//  * EventAction fills them by column NAME (Out::Row<>::set), so there are no numeric column
//    indices anywhere else, and a column cannot be renumbered by accident.
//  * Column names are the branch names in the ROOT file, exactly as written here. tree_chroma and
//    tree_fuse must keep theirs: the GPU reader and Clustering.py look them up by name.
//
//  To add a column: add one X(...) line to the ntuple's list (at the end, to keep old files
//  comparable), then fill it with row.set(Out::<Ntuple>::<name>, value) in EventAction.cc.
//  Types: Int, Dbl, Str.  The ntuple ids below must match the creation order in RunAction.cc
//  (checked at run time).  Units: mm, ns, keV unless the description says otherwise.
// ---------------------------------------------------------------------------------------------
#include "G4AnalysisManager.hh"
#include "G4Exception.hh"
#include "globals.hh"

#include <string>

namespace Out {

// Bump when a column is added, removed or changes meaning; written to run_info as "schema_version".
inline constexpr int kSchemaVersion = 2;

enum class T { Int, Dbl, Str };
struct ColumnSpec { const char* name; T type; const char* description; };

#define OUT_ENUM(type, name, desc) name,
#define OUT_SPEC(type, name, desc) {#name, ::Out::T::type, desc},
#define OUT_NTUPLE(Struct, ID, NAME, TITLE, COLUMNS)                                   \
  struct Struct {                                                                       \
    static constexpr int id = ID;                                                       \
    static constexpr const char* name = NAME;                                           \
    static constexpr const char* title = TITLE;                                         \
    enum Col : int { COLUMNS(OUT_ENUM) nCols };                                         \
    static constexpr ColumnSpec cols[] = { COLUMNS(OUT_SPEC) };                         \
  };

// ---- 0 steps: every step in the LXe (edep > 0 only with /lxe/onlyEdep true) -------------------
#define OUT_COLS_STEPS(X)                                                              \
  X(Int, eventID,    "event number")                                                    \
  X(Int, trackID,    "track id (1 = beam particle)")                                    \
  X(Int, parentID,   "parent track id (0 = primary)")                                   \
  X(Int, pdg,        "PDG code")                                                        \
  X(Dbl, preX_mm,    "step start x")                                                    \
  X(Dbl, preY_mm,    "step start y")                                                    \
  X(Dbl, preZ_mm,    "step start z")                                                    \
  X(Dbl, postX_mm,   "step end x")                                                      \
  X(Dbl, postY_mm,   "step end y")                                                      \
  X(Dbl, postZ_mm,   "step end z")                                                      \
  X(Dbl, edep_keV,   "energy deposit of the step (includes niel)")                      \
  X(Dbl, kinE_keV,   "kinetic energy at step start")                                    \
  X(Dbl, stepLen_mm, "step length")                                                     \
  X(Dbl, time_ns,    "global time at step start")                                       \
  X(Dbl, niel_keV,   "non-ionising part of edep_keV (nuclear stopping)")
OUT_NTUPLE(Steps, 0, "steps", "per-step path points and deposits in LXe", OUT_COLS_STEPS)

// ---- 1 events: one row per event ---------------------------------------------------------------
#define OUT_COLS_EVENTS(X)                                                             \
  X(Int, eventID,        "event number")                                                \
  X(Dbl, edep_keV,       "total deposit in the bulk LXe")                               \
  X(Int, nSteps,         "number of LXe steps stored")                                  \
  X(Dbl, edepPlanes_keV, "deposit inside the fluence planes (not xenon signal)")
OUT_NTUPLE(Events, 1, "events", "per-event totals", OUT_COLS_EVENTS)

// ---- 2 fluence: charged-particle steps in the virtual SiPM planes -----------------------------
#define OUT_COLS_FLUENCE(X)                                                            \
  X(Int, eventID,    "event number")                                                    \
  X(Int, planeID,    "index into Geo::planes / geometry_manifest.json")                 \
  X(Int, trackID,    "track id")                                                        \
  X(Int, parentID,   "parent track id")                                                 \
  X(Int, pdg,        "PDG code")                                                        \
  X(Int, entering,   "1 = step starts on the plane boundary")                           \
  X(Dbl, u_mm,       "local in-plane coordinate u")                                     \
  X(Dbl, v_mm,       "local in-plane coordinate v")                                     \
  X(Dbl, kinE_keV,   "kinetic energy at step start")                                    \
  X(Dbl, cosTheta,   "direction . plane normal")                                        \
  X(Dbl, stepLen_mm, "step length")                                                     \
  X(Dbl, edep_keV,   "deposit in the plane")                                            \
  X(Dbl, time_ns,    "global time")                                                     \
  X(Dbl, vtxX_mm,    "track vertex x (where the particle was created)")                 \
  X(Dbl, vtxY_mm,    "track vertex y")                                                  \
  X(Dbl, vtxZ_mm,    "track vertex z")                                                  \
  X(Str, creator,    "creator process, \"primary\" for the beam particle")
OUT_NTUPLE(Fluence, 2, "fluence", "charged-particle steps in the virtual SiPM planes", OUT_COLS_FLUENCE)

// ---- 3 budget: geantino material budget, one row per (event, volume) --------------------------
#define OUT_COLS_BUDGET(X)                                                             \
  X(Int, eventID,        "event number")                                                \
  X(Str, volume,         "logical volume name")                                         \
  X(Dbl, y0_mm,          "geantino start y")                                            \
  X(Dbl, z0_mm,          "geantino start z")                                            \
  X(Dbl, length_mm,      "path length in the volume")                                   \
  X(Dbl, X0,             "path length / radiation length")                              \
  X(Dbl, lambdaI,        "path length / nuclear interaction length")                    \
  X(Dbl, massThick_gcm2, "mass thickness [g/cm2]")
OUT_NTUPLE(Budget, 3, "budget", "geantino material budget per logical volume", OUT_COLS_BUDGET)

// ---- 4 tree_chroma: input of the GPU optical-photon simulation --------------------------------
// Same layout as the tree_chroma of the other Geant4 simulation that feeds it: one row per
// energy-depositing LXe step, every column a double. "hits" = row index within the event.
// positions and time are post-step; energy and momentum are pre-step. Momentum in keV (p/c).
#define OUT_COLS_CHROMA(X)                                                             \
  X(Dbl, event,         "event number")                                                 \
  X(Dbl, hits,          "row index within the event, 0-based")                          \
  X(Dbl, xpos,          "post-step x")                                                  \
  X(Dbl, ypos,          "post-step y")                                                  \
  X(Dbl, zpos,          "post-step z")                                                  \
  X(Dbl, time,          "post-step global time")                                        \
  X(Dbl, edep,          "energy deposit [keV]")                                         \
  X(Dbl, energy,        "pre-step kinetic energy [keV]")                                \
  X(Dbl, x_momentum,    "pre-step momentum x [keV]")                                    \
  X(Dbl, y_momentum,    "pre-step momentum y [keV]")                                    \
  X(Dbl, z_momentum,    "pre-step momentum z [keV]")                                    \
  X(Dbl, particle_name, "class: 1 gamma, 2 neutron, 3 e-, 4 e+, 5 other (see PDG)")     \
  X(Dbl, PDG,           "PDG code")                                                     \
  X(Dbl, trackID,       "track id")
OUT_NTUPLE(Chroma, 4, "tree_chroma", "Scintillation each hit info", OUT_COLS_CHROMA)

// ---- 5 tree_fuse: the columns Clustering.load_interactions reads ------------------------------
// One row per energy-depositing LXe step plus the zero-deposit gamma interaction steps (they are
// the lineage parents of the electrons they release).
// x, y, z, x_pri, y_pri, z_pri repeat xp.. / xp_pri..: load_interactions asks uproot for both
// spellings in one call and uproot raises on a missing branch.
#define OUT_COLS_FUSE(X)                                                               \
  X(Int, eventid,    "event number")                                                    \
  X(Int, trackid,    "track id")                                                        \
  X(Int, parentid,   "parent track id (0 = primary)")                                   \
  X(Str, type,       "particle name (\"e-\", \"gamma\", \"alpha\", \"Pb208[...]\")")    \
  X(Str, parenttype, "parent's particle name, \"none\" for primaries")                  \
  X(Str, creaproc,   "creator process, \"none\" for primaries; decay -> RadioactiveDecayBase") \
  X(Str, edproc,     "process that ended the step")                                     \
  X(Dbl, ed,         "energy deposit [keV]")                                            \
  X(Dbl, t,          "post-step global time")                                           \
  X(Dbl, xp,         "post-step x")                                                     \
  X(Dbl, yp,         "post-step y")                                                     \
  X(Dbl, zp,         "post-step z")                                                     \
  X(Dbl, xp_pri,     "primary vertex x")                                                \
  X(Dbl, yp_pri,     "primary vertex y")                                                \
  X(Dbl, zp_pri,     "primary vertex z")                                                \
  X(Dbl, x,          "= xp")                                                            \
  X(Dbl, y,          "= yp")                                                            \
  X(Dbl, z,          "= zp")                                                            \
  X(Dbl, x_pri,      "= xp_pri")                                                        \
  X(Dbl, y_pri,      "= yp_pri")                                                        \
  X(Dbl, z_pri,      "= zp_pri")
OUT_NTUPLE(Fuse, 5, "tree_fuse", "per-step lineage table for Clustering.py", OUT_COLS_FUSE)

// ---- 6 primary_path: the beam particle's steps in EVERY volume --------------------------------
#define OUT_COLS_PATH(X)                                                               \
  X(Int, eventID,  "event number")                                                      \
  X(Int, step,     "0 = start point, then the end of each step")                        \
  X(Dbl, x_mm,     "x")                                                                 \
  X(Dbl, y_mm,     "y")                                                                 \
  X(Dbl, z_mm,     "z")                                                                 \
  X(Dbl, kinE_keV, "kinetic energy")                                                    \
  X(Dbl, time_ns,  "global time")                                                       \
  X(Str, volume,   "physical volume the step ran in")
OUT_NTUPLE(Path, 6, "primary_path", "beam particle path through all volumes", OUT_COLS_PATH)

// ---- 7 run_info: provenance of the run, one (key, value) row per item; see RunMeta.cc ----------
#define OUT_COLS_META(X)                                                               \
  X(Str, key,   "item name (cmd_NNN = the NNNth command executed, /vis/ commands left out)") \
  X(Str, value, "item value, as text")
OUT_NTUPLE(Meta, 7, "run_info", "provenance of the run: versions, physics, cuts, seeds, commands", OUT_COLS_META)

#undef OUT_ENUM
#undef OUT_SPEC
#undef OUT_NTUPLE

// ---------------------------------------------------------------------------------------------
//  Creating and filling
// ---------------------------------------------------------------------------------------------
// Creates ntuple N with all its columns, in table order. Call in the order of the ids.
template <class N>
void Create(G4AnalysisManager* man)
{
  const G4int id = man->CreateNtuple(N::name, N::title);
  for (const auto& c : N::cols) {
    switch (c.type) {
      case T::Int: man->CreateNtupleIColumn(c.name); break;
      case T::Dbl: man->CreateNtupleDColumn(c.name); break;
      case T::Str: man->CreateNtupleSColumn(c.name); break;
    }
  }
  man->FinishNtuple();
  if (id != N::id)
    G4Exception("Out::Create", "Out001", FatalException,
                (std::string("ntuple '") + N::name + "' got id " + std::to_string(id) + ", the schema says "
                 + std::to_string(N::id) + ": create the ntuples in id order").c_str());
}

// One row of ntuple N: set columns by name, then add(). A number goes into an Int or a Dbl column
// as the column's type requires; a string only into a Str column.
template <class N>
class Row
{
 public:
  explicit Row(G4AnalysisManager* man) : fMan(man) {}

  void set(typename N::Col c, G4double v) const
  {
    if (N::cols[c].type == T::Int)      fMan->FillNtupleIColumn(N::id, c, static_cast<G4int>(v));
    else if (N::cols[c].type == T::Dbl) fMan->FillNtupleDColumn(N::id, c, v);
    else                                Wrong(c, "a number");
  }
  void set(typename N::Col c, const std::string& v) const
  {
    if (N::cols[c].type == T::Str) fMan->FillNtupleSColumn(N::id, c, v);
    else                           Wrong(c, "a string");
  }
  void set(typename N::Col c, const char* v) const { set(c, std::string(v)); }
  void add() const { fMan->AddNtupleRow(N::id); }

 private:
  static void Wrong(typename N::Col c, const char* what)
  {
    G4Exception("Out::Row::set", "Out002", FatalException,
                (std::string("column '") + N::cols[c].name + "' of ntuple '" + N::name + "' does not take "
                 + what).c_str());
  }
  G4AnalysisManager* fMan;
};

}  // namespace Out

#endif
