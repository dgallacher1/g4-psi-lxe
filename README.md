# Xenon cryostat: Geant4 simulation

Geant4 (>= 11) application: coaxial steel vessels with an insulating-vacuum jacket,
a DN40CF beam port with an Al window, a free-standing sleeve and a non-standard CF on
the inner vessel closed by a 50 um Ti foil into liquid xenon. Scoring in the LXe, plus
two thin virtual planes that record charged-particle fluence where the SiPM boards sit.
No optical physics.

## Frame (matches the GPU light simulations)
- **+z**: vessel axis. **z = 0 = centre of the beam window.**
- **+x**: beam direction. The port is on the -x side; the beam enters from -x.
- x = y = 0 on the vessel axis; y completes the right-handed frame.
- Port parts are located by **s** = distance from the vessel axis along the port (s = -x).

## Geometry: edit `include/GeometryParameters.hh` only
Every number is wrapped as `GEO_V(placeholder)`. Replace the whole wrapper with the CAD value.

| Build option | Unreplaced `GEO_V` entries |
|---|---|
| `-DGEO_DEMO=ON` (default) | placeholders are used; loud DEMO banner; output is `xenon_sim_DEMO.root` |
| `-DGEO_DEMO=OFF` | become NaN; the app **refuses to build** and lists exactly what is missing |

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

## Build & run
```bash
mkdir build && cd build
cmake -DCMAKE_PREFIX_PATH=/path/to/geant4-install [-DGEO_DEMO=OFF] ..
make -j
./geometryProbe                    # geometry self-check
./xenonSim                         # interactive (macros/vis.mac, y=0 section view)
./xenonSim macros/budget_scan.mac  # geantino material budget
./xenonSim macros/run.mac          # GPS positron beam
./xenonSim macros/beam_file.mac    # G4Beamline phase-space file
```
Needs Geant4 built with `GEANT4_USE_GDML=ON` (alpha-source import). Macros live under `macros/`.

## Beam
- `/beam/mode gps`: GPS beam, defined in `macros/beam.mac` (rot1/rot2 put the beam plane normal to +x).
- `/beam/mode file`: G4Beamline BLTrackFile (`x y z Px Py Pz t PDGid EventID TrackID ParentID Weight`).
  Mapped `sim = (bl.z, bl.x, bl.y) + /beam/offset{X,Y,Z}`. One **shared, mutex-protected** reader,
  so MT workers consume disjoint records; the file is rewound each run. Use `/run/beamOn <n records>`.

## Outputs (`xenon_sim.root`)
| ntuple | content |
|---|---|
| `steps` | every step in the bulk LXe (full polylines, edep, PDG, IDs) |
| `events` | per-event `edep_keV` (bulk LXe), `nSteps`, `edepPlanes_keV` |
| `fluence` | every charged-particle step in the planes: plane id, local (u, v), entry KE, cos(theta) to normal, step length, `entering` flag |
| `budget` | geantino scan: per (event, volume) length, X/X0, L/lambda, mass thickness, start (y0, z0) |

Total energy deposited in xenon = `edep_keV + edepPlanes_keV` (the planes are separate LXe volumes).

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
- No optical physics, no hadronic list (MeV-scale e+/alpha).
