#!/usr/bin/env python3
"""Quick look at the LXe response of a (large) run, streamed from tree_chroma.

usage: python lxe_look.py [xenon_sim.root] [geometry_manifest.json] [n_events_to_draw]
deps : uproot numpy matplotlib      (writes lxe_look.png and, if possible, lxe_tracks.png; no display needed)

lxe_look.png   : histograms filled from `events` and `tree_chroma` (every row is an energy deposit in the LXe).
                 tree_chroma is read in chunks into fixed-bin histograms, so memory does not grow with the
                 number of events; a 100k-event file is fine.
lxe_tracks.png : track views from the `steps` ntuple, only when that ntuple is small enough (<= MAX_STEP_ROWS
                 rows) to load; for a big production run use a short run with /output/steps true for this.

Frame: beam along +x, z = vessel axis.  If the manifest is found, the cavity, the Ti foil and the fluence
planes are drawn on the maps.
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

MAX_STEP_ROWS = 3_000_000                                  # above this the track figure is skipped
CHUNK = "300 MB"                                           # uproot.iterate step size

root_file = sys.argv[1] if len(sys.argv) > 1 else "xenon_sim.root"
manifest_file = sys.argv[2] if len(sys.argv) > 2 else "geometry_manifest.json"
n_draw = int(sys.argv[3]) if len(sys.argv) > 3 else 5

geo = json.load(open(manifest_file)) if os.path.exists(manifest_file) else None

PDG = {11: "e-", -11: "e+", 22: "gamma", 2212: "p", 1000020040: "alpha", 13: "mu-", -13: "mu+"}
name = lambda c: PDG.get(int(c), f"ion Z{(int(c) // 10000) % 1000} A{(int(c) // 10) % 1000}" if c > 10**9 else str(int(c)))


def overlay_geometry(panel, b):
    """Draw cavity outline, Ti foil and fluence planes on an (x, b) projection; b is "y" (top) or "z" (side)."""
    if not geo:
        return
    R, i = geo["inner"]["cavityR"], geo["inner"]
    if b == "y":                                           # top view
        panel.add_patch(plt.Circle((0, 0), R, fill=False, color="0.4"))
    else:                                                  # side view
        panel.plot([-R, R, R, -R, -R], [i["cavityZmin"], i["cavityZmin"], i["cavityZmax"], i["cavityZmax"], i["cavityZmin"]], color="0.4")
    panel.axvline(-geo["tiWindow"]["outerFaceS"], color="orange", lw=0.8)
    for p in geo["planes"]:                                # fluence planes (normals along x or z only)
        cx, cy, cz = p["centre"]
        other = cy if b == "y" else cz
        if abs(p["normal"][0]) > 0.99:                     # normal along x: edge-on, a line at x = cx
            h = p["halfV"] if b == "y" else p["halfU"]
            panel.plot([cx, cx], [other - h, other + h], "r-", lw=2)
        elif abs(p["normal"][2]) > 0.99:                   # normal along z: a box (top view) or a line (side view)
            hu, hv = p["halfU"], (p["halfV"] if b == "y" else 0.0)
            panel.plot([cx - hu, cx + hu, cx + hu, cx - hu, cx - hu],
                       [other - hv, other - hv, other + hv, other + hv, other - hv], "r-", lw=2)


# --- events: per-event totals (small) --------------------------------------------------------
with uproot.open(root_file) as f:
    ev = f["events"].arrays(library="np")
    n_ev = len(ev["eventID"])
    n_rows = f["tree_chroma"].num_entries
    n_steps = f["steps"].num_entries if "steps" in f else 0
edep = ev["edep_keV"]

# --- histogram ranges (fixed before streaming) ---------------------------------------------------
R = geo["inner"]["cavityR"] if geo else 150.0
zlo = geo["inner"]["cavityZmin"] if geo else -150.0
zhi = geo["inner"]["cavityZmax"] if geo else 150.0
xe = np.linspace(-R, R, 161)                               # x edges (cavity diameter)
ye = np.linspace(-R, R, 161)
ze = np.linspace(zlo, zhi, 161)
re = np.linspace(0, R, 61)                                 # distance from the beam axis (y = z = 0)

h_x = np.zeros(len(xe) - 1)
h_r = np.zeros(len(re) - 1)
h_xy = np.zeros((len(xe) - 1, len(ye) - 1))
h_xz = np.zeros((len(xe) - 1, len(ze) - 1))
by_pdg = {}                                                # pdg -> deposited keV
n_chunks = 0

with uproot.open(root_file) as f:
    for ch in f["tree_chroma"].iterate(["xpos", "ypos", "zpos", "edep", "PDG"], step_size=CHUNK, library="np"):
        x, y, z, w, pdg = ch["xpos"], ch["ypos"], ch["zpos"], ch["edep"], ch["PDG"].astype(np.int64)
        h_x += np.histogram(x, xe, weights=w)[0]
        h_r += np.histogram(np.hypot(y, z), re, weights=w)[0]
        h_xy += np.histogram2d(x, y, [xe, ye], weights=w)[0]
        h_xz += np.histogram2d(x, z, [xe, ze], weights=w)[0]
        for c in np.unique(pdg):
            by_pdg[int(c)] = by_pdg.get(int(c), 0.0) + w[pdg == c].sum()
        n_chunks += 1

tot = sum(by_pdg.values())
print(f"events {n_ev} | chroma rows {n_rows} ({n_rows / max(n_ev, 1):.0f} per event) | read in {n_chunks} chunk(s)")
print(f"with LXe deposit {np.mean(edep > 0):.3f} | per event: mean {edep.mean():.1f} keV, "
      f"median(>0) {np.median(edep[edep > 0]) if np.any(edep > 0) else 0:.1f} keV, max {edep.max():.1f} keV")
print(f"total deposit in tree_chroma: {tot / n_ev:.1f} keV/event (events ntuple: {edep.mean():.1f}); "
      f"inside the histogram range: {h_x.sum() / n_ev:.1f}")
for c, e in sorted(by_pdg.items(), key=lambda kv: -kv[1])[:6]:
    print(f"  {name(c):>6}: {100 * e / tot:5.1f} % of the deposited energy ({e / n_ev:.1f} keV/event)")
if "edepPlanes_keV" in ev:
    print(f"energy deposited in the silicon planes (not xenon signal): {ev['edepPlanes_keV'].mean():.2f} keV/event")

# --- figure 1: LXe response ---------------------------------------------------------------------
fig, ax = plt.subplots(2, 3, figsize=(16, 9))
xc, rc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (re[1:] + re[:-1])

ax[0, 0].hist(edep[edep > 0], bins=120, log=True)
ax[0, 0].set(xlabel="energy deposited in LXe per event [keV]", ylabel="events", title=f"deposited energy, {n_ev} events")

ax[0, 1].step(xc, h_x / n_ev, where="mid")
ax[0, 1].set(xlabel="x [mm] (beam along +x)", ylabel="keV per event per bin", title="deposit profile along the beam")
if geo:
    xf = -geo["tiWindow"]["outerFaceS"]
    ax[0, 1].axvline(xf, color="k", ls="--", lw=0.8)
    ax[0, 1].text(xf, ax[0, 1].get_ylim()[1], " Ti foil", va="top", fontsize=8)

species = sorted(by_pdg, key=lambda c: -by_pdg[c])[:6]
ax[0, 2].bar([name(c) for c in species], [by_pdg[c] / n_ev for c in species])
ax[0, 2].set(ylabel="keV per event", title="deposit by particle", yscale="log")
ax[0, 2].tick_params(axis="x", labelrotation=45)

for panel, h, ee, b in ((ax[1, 0], h_xy, ye, "y"), (ax[1, 1], h_xz, ze, "z")):
    m = panel.pcolormesh(xe, ee, np.where(h > 0, h / n_ev, np.nan).T, norm=matplotlib.colors.LogNorm(), cmap="viridis")
    fig.colorbar(m, ax=panel, label="keV per event per bin")
    panel.set(xlabel="x [mm]", ylabel=f"{b} [mm]", title=f"deposit map, {'top' if b == 'y' else 'side'} view")
    panel.set_aspect("equal")
    overlay_geometry(panel, b)
    panel.set_xlim(xe[0], xe[-1])
    panel.set_ylim(ee[0], ee[-1])

ax[1, 2].step(rc, h_r / n_ev, where="mid")
ax[1, 2].set(xlabel="distance from the beam axis [mm]", ylabel="keV per event per bin", yscale="log",
             title="transverse spread of the deposit")

fig.tight_layout()
fig.savefig("lxe_look.png", dpi=130)
print("wrote lxe_look.png")

# --- figure 2: tracks, only from a small `steps` ntuple ---------------------------------------------
if n_steps == 0:
    print("no `steps` ntuple entries: no track figure (run with /output/steps true for that)")
elif n_steps > MAX_STEP_ROWS:
    print(f"`steps` has {n_steps} rows (> {MAX_STEP_ROWS}): track figure skipped; use a short run (~200 events) for it")
else:
    with uproot.open(root_file) as f:
        st = f["steps"].arrays(library="np")
    if np.all(st["edep_keV"] > 0):
        print("note: every stored step has edep > 0 (/lxe/onlyEdep true?): trajectories in lxe_tracks.png are not complete")
    shown = [e for e in ev["eventID"] if np.any(st["eventID"] == e)][:n_draw]
    colour = lambda c: {-11: "tab:red", 11: "tab:blue"}.get(int(c), "0.6")
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.5))
    for panel, (a, b, lab) in zip(ax, (("X", "Y", "y"), ("X", "Z", "z"))):
        for e in shown:
            m = st["eventID"] == e
            seg = np.stack([np.stack([st[f"pre{a}_mm"][m], st[f"pre{b}_mm"][m]], 1),
                            np.stack([st[f"post{a}_mm"][m], st[f"post{b}_mm"][m]], 1)], 1)
            cols = ["k" if t == 1 else colour(p) for t, p in zip(st["trackID"][m], st["pdg"][m])]
            panel.add_collection(LineCollection(seg, colors=cols, linewidths=[1.2 if c == "k" else 0.6 for c in cols]))
        panel.set(xlabel="x [mm]", ylabel=f"{lab} [mm]", title=f"tracks in the LXe, {len(shown)} events ({'top' if lab == 'y' else 'side'} view)")
        panel.set_aspect("equal")
        panel.autoscale()
        overlay_geometry(panel, lab)
    fig.tight_layout()
    fig.savefig("lxe_tracks.png", dpi=130)
    print("wrote lxe_tracks.png  (black = beam particle, red e+, blue e-, grey other)")
