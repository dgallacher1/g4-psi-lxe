# Xenon cryostat: Geant4 simulation

Geant4 (>= 11) application: coaxial steel vessels with an insulating-vacuum jacket,
a DN40CF beam port with an Al window, a free-standing sleeve and a non-standard CF on
the inner vessel closed by a 50 um Ti foil into liquid xenon. Scoring in the LXe, plus
two thin silicon planes (`planeMaterial`, G4_Si) that record charged-particle fluence where the
SiPM boards sit.
No optical physics.

## Frame (matches the GPU light simulations)
- **+z**: vessel axis. **z = 0 = centre of the beam window.**
- **+x**: beam direction. The port is on the -x side; the beam enters from -x.
- x = y = 0 on the vessel axis; y completes the right-handed frame.
- Port parts are located by **s** = distance from the vessel axis along the port (s = -x).

## Geometry: edit `include/GeometryParameters.hh` only
All numbers are CAD values now. A new, not yet known dimension is entered as `GEO_V(placeholder)`;
replace the whole wrapper with the CAD value once known.

| Build option | Unreplaced `GEO_V` entries |
|---|---|
| `-DGEO_DEMO=OFF` (**use this**) | become NaN; the app **refuses to build** and lists exactly what is missing. Output is `xenon_sim[_<tag>].root` |
| `-DGEO_DEMO=ON` (CMake default) | placeholders are used; loud DEMO banner; output is `xenon_sim_DEMO[_<tag>].root` |

`GEO_DEMO` defaults to ON in `CMakeLists.txt`, so pass `-DGEO_DEMO=OFF` on the first configure (it is cached in
`CMakeCache.txt`). With no placeholder left in the header, ON only adds the DEMO banner and the `_DEMO` file name.

Each run also writes `geometry_manifest.json` (all dimensions, plane centres/normals/areas,
`isDemo`); the analysis scripts read it, so plane areas can never drift from the C++.

Construction (all solids built in world coordinates, `include/GeometryBuilders.hh`):
- vessel wall = tube minus the port bore; stub = horizontal tube whose root is cut by the
  vertical vessel's outer cylinder; CF flanges are additive rings;
- insulating vacuum and LXe are unions (cavity + port bore) so they reach their foils exactly;
- the Ti foil sits in the inner bore; the sleeve is deliberately unattached.

**After every edit of the header, run `./geometryProbe`.** It validates consistency (flange
clearance, nesting, foils inside their bores) and probes the real solids (port on -x only, bores
cut, stub roots cut, on-axis foil thickness, LXe/vacuum extents, planes inside the LXe).

## Materials
NIST materials, except two custom ones defined in `src/Materials.cc` (the only place to edit a composition or density):
`SS304` (all vessel walls, bottoms, tops, stubs and flanges; nominal 304 by mass: Cr 19, Ni 9.25, Mn 1, Si 0.5, trace C/P/S,
Fe balance, 7.93 g/cm3) and `TiWindowAlloy` (**placeholder**: pure Ti at the G4_Ti density until the alloy is confirmed; the
code prints a PLACEHOLDER line at start-up and flags it in `run_info`). `Materials::Define()` runs in `xenonSim.cc` before the
geometry is built; `GeometryParameters.hh` refers to the names. The `run_info` ntuple records every material used by the
geometry with density and composition (`material[...]`), and which logical volume uses which (`lv[...]`).
The alpha-source GDML still uses `G4_STAINLESS-STEEL` for its puck.

## Build & run
```bash
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=/path/to/geant4-install -DGEO_DEMO=OFF ..
make -j                            # the git commit is regenerated at every make; new macros/*.mac are staged at configure time: re-run cmake after adding one
./geometryProbe                    # geometry self-check
./xenonSim macros/vis.mac          # geometry images (RayTracer, no display needed)
./xenonSim macros/budget_scan.mac  # geantino material budget
./xenonSim macros/smoke_test.mac   # 200 events, checks the output trees (then: python scripts/smoke_check.py xenon_sim_smoke.root)
./xenonSim macros/run.mac          # GPS positron beam
./xenonSim macros/beam_file.mac    # G4Beamline phase-space file
```
Needs Geant4 built with `GEANT4_USE_GDML=ON` (alpha-source import). Macros live under `macros/`.

## Geometry images over SSH (no OpenGL/Qt)
`./xenonSim macros/vis.mac` uses Geant4's built-in RayTracer and writes
`g4RayTracer.RayTracer_NNNN.jpeg` in the working directory (0000 overview, 0001 beam's eye,
0002 outer vessel hidden, 0003 plan view with the top plate hidden). Copy them off with `scp`.
`G4UI_USE_TCSH=1 ./xenonSim` gives a terminal session if you want to issue more `/vis/` commands.
The RayTracer renders geometry only (no trajectories) and ignores section planes, so internals are
exposed with `/vis/geometry/set/visibility <LV> 0 false`. Default image is 640x640; see the
commented `/vis/rayTracer/column|row` lines in the macro.

## Beam
- `/beam/mode gps`: GPS beam, defined in `macros/beam.mac` (rot1/rot2 put the beam plane normal to +x).
- `/beam/mode file`: G4Beamline BLTrackFile (`x y z Px Py Pz t PDGid EventID TrackID ParentID Weight`).
  Mapped `sim = (bl.z, bl.x, bl.y) + /beam/offset{X,Y,Z}`. One **shared, mutex-protected** reader,
  so MT workers consume disjoint records; the file is rewound each run. Use `/run/beamOn <n records>`.

## Outputs (`xenon_sim.root`)
**`include/OutputSchema.hh` is the one place that defines every ntuple and column** (name, type, description).
`RunAction` creates the ntuples from it and `EventAction` fills them by column name (`row.set(Out::Fuse::ed, ...)`),
so there are no numeric column indices in the code. To add a column: one `X(...)` line at the end of the ntuple's list
there, then one `row.set(...)` in the matching `Fill...Row` function of `EventAction.cc`.
Switches (master commands, between runs): `/output/steps`, `/output/chroma`, `/output/fuse`, `/output/path`, `/lxe/onlyEdep`.

| ntuple | content |
|---|---|
| `steps` | every step in the bulk LXe (full polylines, edep, `niel_keV` = non-ionising part of edep, PDG, IDs). `/lxe/onlyEdep true` keeps only edep > 0 rows (smaller; not for trajectory plots) |
| `tree_chroma` | input of the GPU optical-photon simulation, same layout as the other simulation's: one row per energy-depositing LXe step, 14 double columns `event hits xpos ypos zpos time edep energy x_momentum y_momentum z_momentum particle_name PDG trackID` (mm, ns, keV; `hits` = row index within the event; `particle_name` = 1 gamma, 2 neutron, 3 e-, 4 e+, 5 other; positions and time post-step, `energy`/momentum pre-step). `/output/chroma false` and `/output/steps false` switch the trees off |
| `run_info` | provenance of the run, one (key, value) row per item: schema version, UTC time, Geant4 version, compiler, git commit (`-dirty` = uncommitted changes compiled in), geometry fingerprint, events, threads, master seeds, output switches, physics constructors, EM flags, production cuts per region, and every command executed so far (`cmd_NNN`, `/vis/` left out: beam, cuts, seeds, `/lxe/onlyEdep`). The same content is written at the start of each run to `<output>_run_info.json`, so a crashed run still has it |
| `tree_fuse` | lineage table for `Clustering.py`, see below |
| `primary_path` | the beam particle's steps in every volume, see below |
| `events` | per-event `edep_keV` (xenon), `nSteps`, `edepPlanes_keV` (silicon planes) |
| `fluence` | every charged-particle step in the planes: plane id, local (u, v), entry KE, cos(theta) to normal, step length, `entering` flag |
| `budget` | geantino scan: per (event, volume) length, X/X0, L/lambda, mass thickness, start (y0, z0) |

The planes are real silicon, so they scatter and absorb: a plane upstream of another shadows it.

Fluence per primary, per plane (`scripts/fluence.py`):
crossing estimator `sum_entering 1/|cos| / A`, and track-length `sum(stepLen) / (A t)`.

## Analysis
```bash
python scripts/fluence.py      xenon_sim.root geometry_manifest.json
python scripts/budget_map.py   xenon_sim.root geometry_manifest.json
python scripts/quick_look.py   xenon_sim.root
```

## Scope and known limits
- Foil regions use a 2 um production cut (`WindowRegion`). EM option4 handles msc in thin foils;
  validate energy loss through the Ti against expectation for your beam energy.
- No bolt holes, knife edges, gaskets, MLI or cabling; no gas phase above the liquid level.
- Inner CF: no stub; the flange fills a wall opening of radius `flangeOR` and has a cavity-side pocket
  (`innerPocketR`) holding the SiPM, the Ti foil welded at the pocket floor. The LXe is built through the
  opening minus the flange metal. The SiPM holder, PCB, gasket and retaining ring are not modelled.
- Not modelled (off the beam path, small solid angle): the side hole in the inner wall and the lower conflat flange.
- No sleeve in the CAD (`sleeve.use = false`).
- Circular fluence planes (`circular = true`) are discs; the manifest carries each plane's `area`.
- No optical physics, no hadronic list (MeV-scale e+/alpha).

## GPU-simulation input (`tree_chroma`)

`xenon_sim.root` can be fed to `ReadROOTFile(...)` as is. Rows of a multi-threaded run are not in event order
and an event's rows are not guaranteed to be contiguous, so group or sort by `event` (then `hits`).

## Lineage clustering input (`tree_fuse`)

Read by the existing `Clustering.py`; no clustering is done on the G4 side. One row per energy-depositing LXe
step plus the zero-deposit gamma interaction steps (parents of the electrons they release).
Columns: `eventid trackid parentid type parenttype creaproc edproc ed[keV] t[ns] xp yp zp xp_pri yp_pri zp_pri [mm]`
and the same coordinates again as `x y z x_pri y_pri z_pri` (uproot raises on a missing branch, and the loader only keeps
`x_pri` when the plain-named column exists). `creaproc` of any radioactive-decay process is written as `RadioactiveDecayBase`.

    lin, merged = cluster_and_merge_from_root("xenon_sim.root", tree="tree_fuse")

`tree=` is required: with `tree=None` the loader takes the first TTree in the file, which is `steps`.
Switch off with `/output/fuse false`; `/lxe/onlyEdep true` keeps these rows (it only drops zero-deposit steps that are not gamma interactions).

## Beam paths in the real geometry (`primary_path`, `geometryProbe --map`, `beam_path.py`)
No graphics driver or viewer is needed.
- `./geometryProbe --map [step_mm]` samples the solids the application places (`GeoBuild::Build()`, i.e. the CAD numbers of
  `GeometryParameters.hh`) on a 0.5 mm grid in the z = 0 and y = 0 planes and writes `geometry_map_xy.txt` and
  `geometry_map_xz.txt`; it also counts sample points inside two solids that should be disjoint (overlap check, expect 0).
  The alpha-source GDML is not part of it.
- ntuple `primary_path` (switch `/output/path`, default on): the start point and the end of every step of the beam particle
  (trackID 1) in **every** volume (air, Al foil, vacuum, Ti foil, LXe), with the volume name.
- `python scripts/beam_path.py xenon_sim.root geometry_manifest.json [n_events]` draws top and side views of the primaries over
  the map (run it in the directory that holds the two map files; without them it draws from the manifest numbers; without
  `primary_path` it falls back to a dashed straight line from the vertex to the first LXe step).
