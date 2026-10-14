#ifndef RunMeta_hh
#define RunMeta_hh 1

#include <string>
#include <utility>
#include <vector>

class G4Run;

// Provenance of a run: what produced this output file. Collected once on the master at the start of
// each run, written to <output>_run_info.json at once (so a crashed run still has it) and to the
// "run_info" ntuple of the ROOT file at the end of the run.
namespace RunMeta {

using Items = std::vector<std::pair<std::string, std::string>>;

// Master thread, in BeginOfRunAction. `stem` = output file name without ".root".
void Collect(const G4Run* run, const std::string& stem);

// The items collected for the current run (read by the thread that fills the ntuple).
const Items& Get();

}  // namespace RunMeta

#endif
