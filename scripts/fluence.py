#!/usr/bin/env python3
"""Charged-particle fluence at the virtual SiPM planes.

usage: python fluence.py [xenon_sim.root] [geometry_manifest.json]
deps : uproot numpy matplotlib

Per plane, per primary:
  * crossing fluence  Phi_x = sum_{entering rows} 1/|cos(theta)| / A / N_primaries
  * track-length      Phi_l = sum(stepLen) / (A * thickness) / N_primaries
    (they agree for a thin plane; Phi_l also counts particles born in the plane)
plus the entry-energy spectrum, species mix and the (u, v) hit map.
A, thickness come from geometry_manifest.json, written by the sim at construction.
"""
import json
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt

root_file = sys.argv[1] if len(sys.argv) > 1 else "xenon_sim.root"
manifest_file = sys.argv[2] if len(sys.argv) > 2 else "geometry_manifest.json"

with open(manifest_file) as f:
    geo = json.load(f)
if geo["isDemo"]:
    print("*** WARNING: DEMO geometry (placeholder numbers, not CAD) ***\n")

with uproot.open(root_file) as f:
    fl = f["fluence"].arrays(library="np")
    n_primaries = f["events"].num_entries

PDG = {-11: "e+", 11: "e-", 2212: "p", 1000020040: "alpha", 13: "mu-", -13: "mu+"}

fig, axes = plt.subplots(len(geo["planes"]), 3, figsize=(13, 4 * len(geo["planes"])), squeeze=False)
for row, pl in enumerate(geo["planes"]):
    A = 4.0 * pl["halfU"] * pl["halfV"]                  # mm^2
    t = pl["thickness"]
    m = fl["planeID"] == pl["id"]
    ent = m & (fl["entering"] == 1)
    cos = np.abs(fl["cosTheta"][ent])
    ok = cos > 1e-3                                       # guard grazing-angle divergence

    phi_x = np.sum(1.0 / cos[ok]) / A / n_primaries
    phi_l = np.sum(fl["stepLen_mm"][m]) / (A * t) / n_primaries
    n_cross = int(ent.sum())
    print(f"{pl['name']:>6}: area {A:8.1f} mm^2 | crossings {n_cross:8d} "
          f"({n_cross / n_primaries:.4g}/primary) | "
          f"fluence  crossing {phi_x:.4e}  track-length {phi_l:.4e}  [/mm^2/primary]")
    if n_cross:
        stat = np.sqrt(n_cross) / n_cross
        print(f"        stat. rel. error ~ {stat:.2%};  mean E_in = {fl['kinE_keV'][ent].mean():.1f} keV")
        for code in np.unique(fl["pdg"][ent]):
            print(f"        {PDG.get(int(code), int(code)):>6}: {(fl['pdg'][ent] == code).sum()}")

    ax = axes[row]
    ax[0].hist(fl["kinE_keV"][ent], bins=100)
    ax[0].set(xlabel="entry kinetic energy [keV]", ylabel="crossings", title=f"{pl['name']} spectrum")
    ax[1].hist2d(fl["u_mm"][ent], fl["v_mm"][ent], bins=40,
                 range=[[-pl["halfU"], pl["halfU"]], [-pl["halfV"], pl["halfV"]]])
    ax[1].set(xlabel="u [mm]", ylabel="v [mm]", title=f"{pl['name']} hit map", aspect="equal")
    ax[2].hist(cos, bins=50, range=(0, 1))
    ax[2].set(xlabel="|cos(theta)| to plane normal", ylabel="crossings", title=f"{pl['name']} angles")

fig.tight_layout()
fig.savefig("fluence.png", dpi=130)
print("\nwrote fluence.png")
