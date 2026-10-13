#!/usr/bin/env python3
"""Minimal post-processing of xenon_sim.root.

deps:  pip install uproot awkward numpy matplotlib
usage: python quick_look.py [xenon_sim.root]

The g4tools ROOT ntuples map straight onto flat arrays, so uproot reads them
without ROOT installed. 'steps' is the per-step table (path + edep), 'events'
is the per-event summary.
"""
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt

fname = sys.argv[1] if len(sys.argv) > 1 else "xenon_sim.root"

with uproot.open(fname) as f:
    steps = f["steps"].arrays(library="np")
    events = f["events"].arrays(library="np")

# --- per-event energy deposited in the LXe -----------------------------------
edep = events["edep_keV"]
print(f"events: {edep.size} | mean E_dep = {edep.mean():.2f} keV "
      f"| fraction with E_dep>0: {(edep > 0).mean():.3f}")

fig, ax = plt.subplots(1, 2, figsize=(11, 4))

ax[0].hist(edep[edep > 0], bins=120)
ax[0].set_xlabel("energy deposited in LXe per event [keV]")
ax[0].set_ylabel("events")
ax[0].set_title("LXe deposited-energy spectrum")

# --- one event's track polyline in the r-z projection ------------------------
# Pick the event with the most steps for a busy display.
ev = events["eventID"][np.argmax(events["nSteps"])]
m = steps["eventID"] == ev
z = steps["preZ_mm"][m]
r = np.hypot(steps["preX_mm"][m], steps["preY_mm"][m])
c = steps["edep_keV"][m]
sc = ax[1].scatter(z, r, c=c, s=6, cmap="viridis")
ax[1].set_xlabel("z [mm]")
ax[1].set_ylabel("r [mm]")
ax[1].set_title(f"step points in LXe, event {ev}")
fig.colorbar(sc, ax=ax[1], label="step edep [keV]")

fig.tight_layout()
out = "quick_look.png"
fig.savefig(out, dpi=130)
print(f"wrote {out}")
