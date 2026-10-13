#!/usr/bin/env python3
"""Geantino material-budget map: where is the metal in the beam path?

usage: python budget_map.py [xenon_sim.root] [geometry_manifest.json]
deps : uproot numpy matplotlib

Reads the 'budget' ntuple (one row per event x volume, up to first LXe entry),
prints the on-axis and window-aperture-averaged budget per volume, and plots
total X/X0 as a function of the geantino start point (y0, z0).  Compare the
printed per-volume numbers against a hand calculation from the CAD.
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
    b = f["budget"].arrays(library="np")
if b["eventID"].size == 0:
    sys.exit("budget ntuple is empty: run macros/budget_scan.mac first")

vols = np.array([str(v) for v in b["volume"]])
n_events = np.unique(b["eventID"]).size
rap = geo["innerPort"]["boreR"]                       # Ti-window aperture radius
r0 = np.hypot(b["y0_mm"], b["z0_mm"])

print(f"{n_events} geantinos\n")
print(f"{'volume':<14}{'<X/X0> in aperture':>22}{'<L/lambda> in aperture':>26}{'<rho*L> [g/cm2]':>18}")
for v in sorted(set(vols)):
    sel = (vols == v) & (r0 < rap)
    n_in = np.unique(b["eventID"][r0 < rap]).size or 1
    print(f"{v:<14}{b['X0'][sel].sum() / n_in:>22.4e}{b['lambdaI'][sel].sum() / n_in:>26.4e}"
          f"{b['massThick_gcm2'][sel].sum() / n_in:>18.4e}")

# ---- total X/X0 vs start point ------------------------------------------------
ev = np.unique(b["eventID"])
tot = np.zeros(ev.size)
y0 = np.zeros(ev.size)
z0 = np.zeros(ev.size)
idx = np.searchsorted(ev, b["eventID"])
np.add.at(tot, idx, b["X0"])
y0[idx] = b["y0_mm"]
z0[idx] = b["z0_mm"]

fig, ax = plt.subplots(figsize=(6.5, 5.5))
h = ax.hexbin(y0, z0, C=tot, reduce_C_function=np.mean, gridsize=60)
fig.colorbar(h, ax=ax, label="mean total X/X0 before LXe")
ax.set(xlabel="start y [mm]", ylabel="start z [mm]", title="Material budget to the LXe", aspect="equal")
fig.tight_layout()
fig.savefig("budget_map.png", dpi=130)
print("\nwrote budget_map.png")
