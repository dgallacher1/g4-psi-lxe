# Xenon cryostat — basic Geant4 simulation

A minimal but complete Geant4 (>= 11) application: a concentric cylindrical
cryostat holding liquid xenon, a GDML/CAD-imported alpha source submerged in
the LXe, a runtime-configurable positron beam, and step-level scoring in the
xenon. No optical physics.

Author - David  G& Claude code

## Requirements
- Geant4 >= 11.0, built with `GEANT4_USE_GDML=ON` (Xerces-C).
- Multithread build recommended (analysis ntuple merging is enabled).
- Vis (OpenGL/Qt) only for the interactive viewer.

## Build & run
```bash
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=/home/dgallacher/packages/geant4/g4-install/v11.4.2 ..
make -j
./xenonSim              # interactive viewer (macros/vis.mac)
./xenonSim run.mac      # batch, writes xenon_sim.root
```
Assets (`macros/`, `gdml/`) are staged into the build dir by CMake, so relative
paths work from there.

## Geometry
`world(vac) > outerVessel(steel) > vacuumJacket(vac) > innerVessel(steel) > LXe`.
Nested *solid* `G4Tubs` — each daughter carves its parent, so the steel vessels
are closed (walls + end caps) without boolean solids. Edit dimensions at the top
of `DetectorConstruction::Construct()`.

The alpha source is read from `gdml/alpha_source.gdml`; the code extracts the LV
named `AlphaSourceLV` and places it at the LXe centre. Swap in your CAD export
(STEP -> GDML mesh); keep the top LV name or change `fSourceLVName`.

## Beam
- `mode gps` (default): positron beam via GPS, fully defined in `macros/beam.mac`
  (spot, direction, spectrum, optional divergence). No recompilation to change it.
- `mode file`: read a real G4Beamline phase-space dump (BLTrackFile:
  `x y z Px Py Pz t PDGid EventID TrackID ParentID Weight`, mm/MeV/ns), one
  primary per event:
  ```
  /beam/mode file
  /beam/file  positrons_at_cryostat.txt
  /run/beamOn 0     # runs until the file is exhausted
  ```

## Scoring / output
`XenonSD` on the LXe LV records every step in the xenon as a hit (pre/post point,
edep, PDG, track/parent ID, KE, step length, time). `EventAction` streams these to
two ntuples in `xenon_sim.root`:
- `steps`  — one row per step (full track polylines + deposits)
- `events` — one row per event (total edep, step multiplicity)

Switch to CSV/HDF5 with `SetDefaultFileType` in `RunAction`. Analyse with
`analysis/quick_look.py` (uproot).

## Deliberate scope limits
- No `G4OpticalPhysics`: scintillation/SiPM response is out of scope; add it plus a
  photon-sensitive SD if you later want light readout.
- No hadronic list: fine for MeV-scale e+/alpha.
