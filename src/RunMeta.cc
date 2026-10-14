#include "RunMeta.hh"
#include "GeometryParameters.hh"
#include "OutputOptions.hh"
#include "OutputSchema.hh"

#include "G4Version.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4EmParameters.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <cstdint>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>

#if __has_include("GitVersion.hh")
#include "GitVersion.hh"                 // generated at build time by cmake/GitVersion.cmake
#endif
#ifndef XENON_GIT_DESCRIBE
#define XENON_GIT_DESCRIBE "unknown"
#endif

namespace RunMeta {
namespace {

Items gItems;      // written by the master at the start of a run, read by workers at its end

std::string Hex64(std::uint64_t v)
{
  std::ostringstream o;
  o << std::hex << std::setw(16) << std::setfill('0') << v;
  return o.str();
}

// FNV-1a: a short fingerprint of the geometry manifest (all dimensions), stable across platforms.
std::uint64_t Fnv1a(const std::string& s)
{
  std::uint64_t h = 1469598103934665603ULL;
  for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
  return h;
}

std::string Utc()
{
  const std::time_t t = std::time(nullptr);
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return buf;
}

std::string JsonEscape(const std::string& s)
{
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\') { o += '\\'; o += c; }
    else if (c == '\n') o += "\\n";
    else if (c == '\t') o += "\\t";
    else o += c;
  }
  return o;
}

std::string Num(double v, int prec = 6)
{
  std::ostringstream o;
  o << std::setprecision(prec) << v;
  return o.str();
}

}  // namespace

void Collect(const G4Run* run, const std::string& stem)
{
  Items it;
  auto add = [&](const std::string& k, const std::string& v) { it.emplace_back(k, v); };
  auto yn = [](bool b) { return std::string(b ? "true" : "false"); };

  // ---- software and build ----
  add("schema_version", std::to_string(Out::kSchemaVersion));
  add("timestamp_utc", Utc());
  add("geant4_version", G4Version);
  add("compiler", __VERSION__);
  add("git", XENON_GIT_DESCRIBE);                 // "-dirty" = uncommitted changes were compiled in
  add("geo_demo", yn(Geo::isDemo));
  add("geometry_fingerprint", Hex64(Fnv1a(Geo::ManifestJSON())));

  // ---- run setup ----
  auto* rm = G4RunManager::GetRunManager();
  add("output_stem", stem);
  add("events_requested", std::to_string(run->GetNumberOfEventToBeProcessed()));
  int nThreads = 1;
  if (auto* mt = dynamic_cast<G4MTRunManager*>(rm)) nThreads = mt->GetNumberOfThreads();
  add("threads", std::to_string(nThreads));
  add("rng_engine", G4Random::getTheEngine()->name());
  if (const long* s = CLHEP::HepRandom::getTheSeeds())
    add("master_seeds", std::to_string(s[0]) + " " + std::to_string(s[1]));
  add("write_steps", yn(gWriteSteps));
  add("write_chroma", yn(gWriteChroma));
  add("write_fuse", yn(gWriteFuse));
  add("write_path", yn(gWritePath));

  // ---- physics ----
  if (auto* pl = dynamic_cast<const G4VModularPhysicsList*>(rm->GetUserPhysicsList())) {
    std::string names;
    for (int i = 0; const auto* c = pl->GetPhysics(i); ++i) names += (i ? ", " : "") + c->GetPhysicsName();
    add("physics_constructors", names);
    add("default_cut_mm", Num(pl->GetDefaultCutValue() / mm));
  }
  auto* em = G4EmParameters::Instance();
  add("em_fluo", yn(em->Fluo()));
  add("em_auger", yn(em->Auger()));
  add("em_pixe", yn(em->Pixe()));
  add("em_deexcitation_ignore_cut", yn(em->DeexcitationIgnoreCut()));
  add("em_lowest_electron_energy_keV", Num(em->LowestElectronEnergy() / keV));
  for (const auto* region : *G4RegionStore::GetInstance()) {
    if (const auto* pc = region->GetProductionCuts()) {
      add("cuts_mm[" + region->GetName() + "]  gamma e- e+ proton",
          Num(pc->GetProductionCut("gamma") / mm) + " " + Num(pc->GetProductionCut("e-") / mm) + " " +
          Num(pc->GetProductionCut("e+") / mm) + " " + Num(pc->GetProductionCut("proton") / mm));
    }
  }

  // ---- every command executed so far (macros included): beam, cuts, seeds, /lxe/onlyEdep ... ----
  auto* ui = G4UImanager::GetUIpointer();
  int n = 0;
  for (int i = 0; i < ui->GetNumberOfHistory(); ++i) {
    const std::string cmd = ui->GetPreviousCommand(i);
    if (cmd.rfind("/vis/", 0) == 0) continue;
    char key[16];
    std::snprintf(key, sizeof key, "cmd_%03d", n++);
    add(key, cmd);
  }

  gItems = it;

  std::ofstream js(stem + "_run_info.json");
  js << "{\n";
  for (std::size_t i = 0; i < it.size(); ++i)
    js << "  \"" << JsonEscape(it[i].first) << "\": \"" << JsonEscape(it[i].second) << "\""
       << (i + 1 < it.size() ? "," : "") << "\n";
  js << "}\n";
}

const Items& Get() { return gItems; }

}  // namespace RunMeta
