#!/usr/bin/env bash
# Run the focused-beam LXe simulation in a fresh directory on /scratch, then lxe_look.py and fluence.py.
#
#   bash scripts/run.sh                      # macros/lxe_focused.mac
#   bash scripts/run.sh macros/other.mac     # another run macro (it must contain /output/tag <tag>)
#   nohup bash scripts/run.sh > run_all.log 2>&1 &     # leave it running; follow with: tail -f run_all.log
#
# Everything the run writes (ROOT file, run_info.json, geometry_manifest.json, PNGs, logs) lands in
#   $SCRATCH_ROOT/<tag>_<date>_<time>/
# Re-run only the analysis for an existing run:   SKIP_SIM=1 OUT=/scratch/.../<dir> bash scripts/run.sh
#
# Settings (environment variables, defaults in brackets):
#   REPO [~/packages/geant4/g4-psi]   BUILD [$REPO/build]   SCRATCH_ROOT [/scratch/dgallacher/psi-beam-sims]
#   G4SH [~/packages/geant4/g4-install/v11.4.2/bin/geant4.sh]   CONDA_ENV [muon_lgad]   MIN_FREE_GB [60]
set -eo pipefail

REPO=${REPO:-$HOME/packages/geant4/g4-psi}
BUILD=${BUILD:-$REPO/build}
SCRATCH_ROOT=${SCRATCH_ROOT:-/scratch/dgallacher/psi-beam-sims}
G4SH=${G4SH:-$HOME/packages/geant4/g4-install/v11.4.2/bin/geant4.sh}
CONDA_ENV=${CONDA_ENV:-muon_lgad}
MIN_FREE_GB=${MIN_FREE_GB:-60}
PY=${PY:-python}
SKIP_SIM=${SKIP_SIM:-0}
MACRO=${1:-macros/lxe_focused.mac}

die() { echo "run.sh: $*" >&2; exit 1; }

[ -f "$REPO/$MACRO" ] || die "macro not found: $REPO/$MACRO"
TAG=$(awk '$1 == "/output/tag" { print $2; exit }' "$REPO/$MACRO")
[ -n "$TAG" ] || die "no '/output/tag' line in $MACRO (needed to know the ROOT file name)"
ROOT=xenon_sim_$TAG.root

if [ "$SKIP_SIM" = 1 ]; then
  [ -n "$OUT" ] && [ -d "$OUT" ] || die "SKIP_SIM=1 needs OUT=<existing run directory>"
else
  OUT=${OUT:-$SCRATCH_ROOT/${TAG}_$(date +%Y%m%d_%H%M%S)}
  [ -x "$BUILD/xenonSim" ]            || die "executable not found: $BUILD/xenonSim"
  [ -f "$REPO/gdml/alpha_source.gdml" ] || die "missing $REPO/gdml/alpha_source.gdml (the sim would run without the source)"
  mkdir -p "$OUT"
  FREE_GB=$(df -Pk "$OUT" | awk 'NR == 2 { printf "%d", $4 / 1048576 }')
  if [ "$FREE_GB" -lt "$MIN_FREE_GB" ]; then
    rmdir "$OUT" 2>/dev/null || true
    die "only ${FREE_GB} GB free on the filesystem of $OUT (MIN_FREE_GB=$MIN_FREE_GB)"
  fi
fi
cd "$OUT"
echo "run directory : $OUT"
echo "macro / file  : $MACRO -> $ROOT"
echo "started       : $(date)"
SECONDS=0

# ---- 1. simulation (before conda is activated, so the Geant4 libraries come from the usual environment) ----
if [ "$SKIP_SIM" != 1 ]; then
  # The sim reads ./macros and ./gdml relative to the working directory and writes its files there.
  ln -sfn "$REPO/macros" macros
  ln -sfn "$REPO/gdml"   gdml
  if [ -f "$G4SH" ]; then . "$G4SH" || true; fi

  "$BUILD/xenonSim" "$MACRO" 2>&1 | tee run.log
  echo "simulation finished after ${SECONDS} s"

  grep -q "running without source" run.log && echo "WARNING: the log says the alpha source was not loaded (check the gdml link)" >&2
  [ -s "$ROOT" ] || die "no $ROOT after the run; see $OUT/run.log"
  ls -lh "$ROOT"
fi
[ -s "$ROOT" ]                    || die "$OUT/$ROOT not found"
[ -f geometry_manifest.json ]     || die "$OUT/geometry_manifest.json not found"

# ---- 2. analysis (needs the conda environment) ----
CONDA_SH=${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}
if [ ! -f "$CONDA_SH" ] && command -v conda >/dev/null 2>&1; then
  CONDA_SH="$(conda info --base)/etc/profile.d/conda.sh"
fi
if [ ! -f "$CONDA_SH" ]; then
  echo "conda.sh not found; the simulation output is in $OUT. Run the analysis by hand:" >&2
  echo "  conda activate $CONDA_ENV && cd $OUT && python $REPO/scripts/lxe_look.py $ROOT geometry_manifest.json && python $REPO/scripts/fluence.py $ROOT geometry_manifest.json" >&2
  exit 1
fi
. "$CONDA_SH"
conda activate "$CONDA_ENV"

echo; echo "=== lxe_look.py ==="
$PY "$REPO/scripts/lxe_look.py" "$ROOT" geometry_manifest.json 2>&1 | tee lxe_look.out
echo; echo "=== fluence.py ==="
$PY "$REPO/scripts/fluence.py" "$ROOT" geometry_manifest.json 2>&1 | tee fluence.out

echo; echo "done after ${SECONDS} s; files in $OUT:"
ls -lh "$OUT"
