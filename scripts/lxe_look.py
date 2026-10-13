#!/usr/bin/env python3
"""Quick look at what the LXe sensitive detector recorded: energy deposits and tracks.

usage: python lxe_look.py [xenon_sim.root] [geometry_manifest.json] [n_events_to_draw]
deps : uproot numpy matplotlib      (writes lxe_look.png; no display needed)

Reads the 'steps' (every step in the LXe) and 'events' ntuples.  Frame: beam along +x,
z = vessel axis.  If the manifest is found, the cavity, the Ti foil and the fluence
planes are drawn on the track views.
"""
import json
import os
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")                                     # works over SSH
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

root_file = sys.argv[1] if len(sys.argv) > 1 else "xenon_sim.root"
manifest_file = sys.argv[2] if len(sys.argv) > 2 else "geometry_manifest.json"
n_draw = int(sys.argv[3]) if len(sys.argv) > 3 else 5

geo = json.load(open(manifest_file)) if os.path.exists(manifest_file) else None
with uproot.open(root_file) as f:
    st = f["steps"].arrays(library="np")
    ev = f["events"].arrays(library="np")
n_ev = len(ev["eventID"])
if len(st["eventID"]) == 0:
    sys.exit("no steps in the LXe in this file (a geantino scan, or the beam missed the LXe?)")

PDG = {11: "e-", -11: "e+", 22: "gamma", 2212: "p", 1000020040: "alpha", 13: "mu-", -13: "mu+"}
name = lambda c: PDG.get(int(c), f"ion Z{(int(c) // 10000) % 1000} A{(int(c) // 10) % 1000}" if c > 10**9 else str(int(c)))

# --- text summary ----------------------------------------------------------------
edep = ev["edep_keV"]
print(f"events {n_ev} | with LXe deposit {np.mean(edep > 0):.3f} | per event: mean {edep.mean():.1f} keV, "
      f"median(>0) {np.median(edep[edep > 0]):.1f} keV, max {edep.max():.1f} keV")
tot = st["edep_keV"].sum()
by_pdg = {int(c): st["edep_keV"][st["pdg"] == c].sum() for c in np.unique(st["pdg"])}
for c, e in sorted(by_pdg.items(), key=lambda kv: -kv[1])[:6]:
    print(f"  {name(c):>6}: {100 * e / tot:5.1f} % of the deposited energy ({e / n_ev:.1f} keV/event)")
if "niel_keV" in st:                                       # non-ionising share (nuclear stopping)
    print(f"non-ionising part of the deposit: {100 * st['niel_keV'].sum() / tot:.4f} %")
if np.all(st["edep_keV"] > 0):
    print("note: every stored step has edep > 0 (/lxe/onlyEdep true?): trajectories below are not complete")
prim = st["trackID"] == 1
n_with_primary = len(np.unique(st["eventID"][prim]))
print(f"primary enters the LXe in {n_with_primary} / {n_ev} events; "
      f"mean primary path in LXe {st['stepLen_mm'][prim].sum() / max(n_with_primary, 1):.1f} mm")

# --- figure ----------------------------------------------------------------------
fig, ax = plt.subplots(2, 3, figsize=(16, 9))

ax[0, 0].hist(edep[edep > 0], bins=120, log=True)
ax[0, 0].set(xlabel="energy deposited in LXe per event [keV]", ylabel="events", title="deposited energy")

xm = 0.5 * (st["preX_mm"] + st["postX_mm"])                # where along the beam the energy goes
ax[0, 1].hist(xm, bins=100, weights=st["edep_keV"] / n_ev)
ax[0, 1].set(xlabel="x [mm] (beam along +x)", ylabel="keV per event per bin", title="deposit profile along the beam")
if geo:
    ax[0, 1].axvline(-geo["tiWindow"]["outerFaceS"], color="k", ls="--", lw=0.8)
    ax[0, 1].text(-geo["tiWindow"]["outerFaceS"], ax[0, 1].get_ylim()[1], " Ti foil", va="top", fontsize=8)

species = sorted(by_pdg, key=lambda c: -by_pdg[c])[:6]
ax[0, 2].bar([name(c) for c in species], [by_pdg[c] / n_ev for c in species])
ax[0, 2].set(ylabel="keV per event", title="deposit by particle", yscale="log")

# tracks of the first few events with LXe steps: black = beam particle, red e+, blue e-, grey other
shown = [e for e in ev["eventID"] if np.any(st["eventID"] == e)][:n_draw]
colour = lambda c: {-11: "tab:red", 11: "tab:blue"}.get(int(c), "0.6")
for panel, (a, b, lab) in zip((ax[1, 0], ax[1, 1]), (("X", "Y", "y"), ("X", "Z", "z"))):
    for e in shown:
        m = st["eventID"] == e
        seg = np.stack([np.stack([st[f"pre{a}_mm"][m], st[f"pre{b}_mm"][m]], 1),
                        np.stack([st[f"post{a}_mm"][m], st[f"post{b}_mm"][m]], 1)], 1)
        cols = ["k" if t == 1 else colour(p) for t, p in zip(st["trackID"][m], st["pdg"][m])]
        panel.add_collection(LineCollection(seg, colors=cols, linewidths=[1.2 if c == "k" else 0.6 for c in cols]))
    panel.set(xlabel="x [mm]", ylabel=f"{lab} [mm]", title=f"tracks in the LXe, {len(shown)} events ({'top' if lab == 'y' else 'side'} view)")
    panel.set_aspect("equal")
    panel.autoscale()
    if geo:
        R, i = geo["inner"]["cavityR"], geo["inner"]
        if lab == "y":
            panel.add_patch(plt.Circle((0, 0), R, fill=False, color="0.4"))
        else:
            panel.plot([-R, R, R, -R, -R], [i["cavityZmin"], i["cavityZmin"], i["cavityZmax"], i["cavityZmax"], i["cavityZmin"]], color="0.4")
        panel.axvline(-geo["tiWindow"]["outerFaceS"], color="orange", lw=0.8)
        for p in geo["planes"]:                             # fluence planes (normals along x or z only)
            cx, cy, cz = p["centre"]
            other = cy if lab == "y" else cz
            if abs(p["normal"][0]) > 0.99:                  # normal along x: edge-on, a line at x = cx
                h = p["halfV"] if lab == "y" else p["halfU"]
                panel.plot([cx, cx], [other - h, other + h], "r-", lw=2)
            elif abs(p["normal"][2]) > 0.99:                # normal along z: a box (top view) or a line (side view)
                hu, hv = p["halfU"], (p["halfV"] if lab == "y" else 0.0)
                panel.plot([cx - hu, cx + hu, cx + hu, cx - hu, cx - hu],
                           [other - hv, other - hv, other + hv, other + hv, other - hv], "r-", lw=2)

r_perp = np.hypot(st["preY_mm"], st["preZ_mm"])            # distance from the beam axis (y = z = 0)
ax[1, 2].hist(r_perp, bins=60, weights=st["edep_keV"] / n_ev, log=True)
ax[1, 2].set(xlabel="distance from the beam axis [mm]", ylabel="keV per event per bin", title="transverse spread of the deposit")

fig.tight_layout()
fig.savefig("lxe_look.png", dpi=130)
print("wrote lxe_look.png")
