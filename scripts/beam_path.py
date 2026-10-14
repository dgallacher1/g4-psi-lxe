#!/usr/bin/env python3
"""Beam paths into the cryostat: top (x-y) and side (x-z) views of the PRIMARY tracks only.

usage: python beam_path.py [xenon_sim.root] [geometry_manifest.json] [n_events=30]
deps : uproot numpy matplotlib      (writes beam_path.png; no display needed)

* Background: the REAL Geant4 solids if geometry_map_xy.txt / geometry_map_xz.txt are in the working
  directory (made by `./geometryProbe --map`, which samples the same GeoBuild::Build() solids the
  application places); otherwise a drawing from the manifest numbers.
* Black lines: the beam particle (trackID 1). With the 'primary_path' ntuple this is its true path in
  every volume (air, Al foil, vacuum, Ti foil, LXe). Without it only the LXe part is known (from
  'steps') and the flight in from the vertex is drawn as a dashed straight line (not tracked).
* Secondaries are not drawn.
"""
import json
import os
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle, Wedge, Circle
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm

root_file = sys.argv[1] if len(sys.argv) > 1 else "xenon_sim.root"
manifest_file = sys.argv[2] if len(sys.argv) > 2 else "geometry_manifest.json"
n_draw = int(sys.argv[3]) if len(sys.argv) > 3 else 30
geo = json.load(open(manifest_file))

STEEL, LXE, VAC = "#b9bec6", "#cfe6f5", "#ffffff"

with uproot.open(root_file) as f:
    st = f["steps"].arrays(["eventID", "trackID", "preX_mm", "preY_mm", "preZ_mm",
                            "postX_mm", "postY_mm", "postZ_mm"], library="np")
    vtx = None
    pp = None
    if "primary_path" in f:
        pp = f["primary_path"].arrays(["eventID", "step", "x_mm", "y_mm", "z_mm"], library="np")
    if "tree_fuse" in f:
        fz = f["tree_fuse"].arrays(["eventid", "xp_pri", "yp_pri", "zp_pri"], library="np")
        ids, first = np.unique(fz["eventid"], return_index=True)
        vtx = {int(e): (fz["xp_pri"][i], fz["yp_pri"][i], fz["zp_pri"][i]) for e, i in zip(ids, first)}
if len(st["eventID"]) == 0:
    sys.exit("no LXe steps in this file")

prim = st["trackID"] == 1
events = [int(e) for e in np.unique(st["eventID"][prim])][:n_draw]
if not events:
    sys.exit("the primary never reaches the LXe in this file")
x0_default = geo.get("suggestedBeamStartX", -175.0)

o, i, op, ip = geo["outer"], geo["inner"], geo["outerPort"], geo["innerPort"]
al, ti = geo["alWindow"], geo["tiWindow"]


LABEL_COL = ["#ffffff", "#f1f1f1", "#cfe6f5", "#b9bec6", "#b9bec6", "#9aa1ab", "#9aa1ab",
             "#9aa1ab", "#2ca02c", "#ff7f0e", "#d9b38c", "#9aa1ab"]


def load_map(fn):
    if not os.path.exists(fn):
        return None
    with open(fn) as fh:
        x0, x1, v0, v1, step, nx, nv = [float(t) for t in fh.readline().split()[1:]]
    return np.loadtxt(fn, dtype=int, skiprows=1), (x0, x1, v0, v1)


maps = {"top": load_map("geometry_map_xy.txt"), "side": load_map("geometry_map_xz.txt")}
use_map = all(m is not None for m in maps.values())


def draw_map(ax, view):
    img, ext = maps[view]
    ax.imshow(img, origin="lower", extent=ext, cmap=ListedColormap(LABEL_COL),
              norm=BoundaryNorm(np.arange(-0.5, len(LABEL_COL)), len(LABEL_COL)), interpolation="nearest", zorder=1)


def port_parts(ax, view):
    """Port stubs, flanges and the inner-vessel port: identical in both views (axisymmetric about x)."""
    def rect(x0, x1, r0, r1, col=STEEL, z=3):          # mirrored about the axis: r in [r0,r1] and [-r1,-r0]
        for s in (1, -1):
            ax.add_patch(Rectangle((x0, min(s * r0, s * r1)), x1 - x0, abs(r1 - r0), fc=col, ec="0.35", lw=0.6, zorder=z))
    # outer vessel: bore cut through the wall (vacuum), stub and flange outside
    ax.add_patch(Rectangle((-(o["cavityR"] + o["wallT"]) - 0.5, -op["boreR"]), o["wallT"] + 1.0, 2 * op["boreR"],
                           fc=VAC, ec="none", zorder=2))
    rect(-op["flangeS0"], -(o["cavityR"] + o["wallT"]), op["boreR"], op["stubOR"])
    rect(-op["flangeS1"], -op["flangeS0"], op["boreR"], op["flangeOR"])
    # inner vessel: LXe-filled bore up to the Ti foil, flange around it
    ax.add_patch(Rectangle((-ti["outerFaceS"], -ip["boreR"]), ti["outerFaceS"] - i["cavityR"] + 0.5, 2 * ip["boreR"],
                           fc=LXE, ec="none", zorder=2.5))
    ax.add_patch(Rectangle((-(i["cavityR"] + i["wallT"]) - 0.01, -ip["boreR"]), i["wallT"] + 0.02, 2 * ip["boreR"],
                           fc=LXE, ec="none", zorder=2.5))
    rect(-ip["flangeS1"], -ip["flangeS0"], ip["boreR"], ip["flangeOR"])
    # foils (thin: drawn as lines)
    ax.plot([-ti["outerFaceS"]] * 2, [-ip["boreR"], ip["boreR"]], color="tab:orange", lw=1.6, zorder=6)
    ax.plot([-al["outerFaceS"]] * 2, [-op["boreR"], op["boreR"]], color="tab:green", lw=1.6, zorder=6)


def draw_top(ax):
    Ro, Ri = o["cavityR"], i["cavityR"]
    ax.add_patch(Circle((0, 0), Ro + o["wallT"], fc=STEEL, ec="0.35", lw=0.6, zorder=1))
    ax.add_patch(Circle((0, 0), Ro, fc=VAC, ec="0.35", lw=0.6, zorder=1.2))
    ax.add_patch(Circle((0, 0), Ri + i["wallT"], fc=STEEL, ec="0.35", lw=0.6, zorder=1.4))
    ax.add_patch(Circle((0, 0), Ri, fc=LXE, ec="0.35", lw=0.6, zorder=1.6))
    port_parts(ax, "top")


def draw_side(ax):
    ax.add_patch(Rectangle((-(o["cavityR"] + o["wallT"]), o["cavityZmin"] - o["bottomT"]),
                           2 * (o["cavityR"] + o["wallT"]), o["cavityZmax"] - o["cavityZmin"] + o["bottomT"] + o["topT"],
                           fc=STEEL, ec="0.35", lw=0.6, zorder=1))
    ax.add_patch(Rectangle((-o["cavityR"], o["cavityZmin"]), 2 * o["cavityR"], o["cavityZmax"] - o["cavityZmin"],
                           fc=VAC, ec="0.35", lw=0.6, zorder=1.2))
    ax.add_patch(Rectangle((-(i["cavityR"] + i["wallT"]), i["cavityZmin"] - i["bottomT"]),
                           2 * (i["cavityR"] + i["wallT"]), i["cavityZmax"] - i["cavityZmin"] + i["bottomT"] + i["topT"],
                           fc=STEEL, ec="0.35", lw=0.6, zorder=1.4))
    ax.add_patch(Rectangle((-i["cavityR"], i["cavityZmin"]), 2 * i["cavityR"], i["cavityZmax"] - i["cavityZmin"],
                           fc=LXE, ec="0.35", lw=0.6, zorder=1.6))
    port_parts(ax, "side")


def planes(ax, view):
    for p in geo["planes"]:
        cx, cy, cz = p["centre"]
        other = cy if view == "top" else cz
        if abs(p["normal"][0]) > 0.99:
            h = p["halfV"] if view == "top" else p["halfU"]
            ax.plot([cx, cx], [other - h, other + h], color="red", lw=2.5, zorder=7)
        elif abs(p["normal"][2]) > 0.99:
            hu, hv = p["halfU"], (p["halfV"] if view == "top" else 0.0)
            ax.plot([cx - hu, cx + hu, cx + hu, cx - hu, cx - hu],
                    [other - hv, other - hv, other + hv, other + hv, other - hv], color="red", lw=2.5, zorder=7)


fig, axs = plt.subplots(2, 1, figsize=(11, 13), sharex=True, gridspec_kw={"height_ratios": [1.0, 1.95]})
for ax, view, (b, lab) in zip(axs, ("top", "side"), (("Y", "y"), ("Z", "z"))):
    if use_map:
        draw_map(ax, view)
    else:
        (draw_top if view == "top" else draw_side)(ax)
    planes(ax, view)
    segs_in, segs_out, ends, starts, entries = [], [], [], [], []
    for e in events:
        m = (st["eventID"] == e) & prim
        pre = np.stack([st["preX_mm"][m], st[f"pre{b}_mm"][m]], 1)
        post = np.stack([st["postX_mm"][m], st[f"post{b}_mm"][m]], 1)
        if len(pre) == 0:
            continue
        entries.append(pre[0])
        ends.append(post[-1])
        mp = (pp["eventID"] == e) if pp is not None else None
        if mp is not None and mp.any():                       # true path through every volume
            order = np.argsort(pp["step"][mp])
            line = np.stack([pp["x_mm"][mp][order], pp[f"{b.lower()}_mm"][mp][order]], 1)
            segs_in.append(np.stack([line[:-1], line[1:]], 1))
            starts.append(line[0])
            continue
        segs_in.append(np.stack([pre, post], 1))
        v = vtx.get(e) if vtx else None
        start = np.array([v[0], v[1] if view == "top" else v[2]]) if v else np.array([x0_default, 0.0])
        starts.append(start)
        segs_out.append(np.stack([start, pre[0]])[None])
    ax.add_collection(LineCollection(np.concatenate(segs_in), colors="k", linewidths=0.9, alpha=0.8, zorder=5))
    if segs_out:
        ax.add_collection(LineCollection(np.concatenate(segs_out), colors="tab:purple", linewidths=0.8,
                                         linestyles="--", alpha=0.6, zorder=4.5))
    s, en, nd = np.array(starts), np.array(entries), np.array(ends)
    ax.scatter(*s.T, marker="*", s=140, c="gold", ec="k", lw=0.7, zorder=8)
    ax.scatter(*en.T, marker="o", s=14, c="tab:purple", zorder=8)
    ax.scatter(*nd.T, marker="x", s=18, c="tab:red", zorder=8)
    ax.set_ylabel(f"{lab} [mm]")
    ax.set_aspect("equal")
    ax.set_title("top view (x-y)" if view == "top" else "side view (x-z)", fontsize=10)
    ax.grid(alpha=0.2)

axs[1].set_xlabel("x [mm]   (beam along +x, port on the -x side)")
xmin = min(-(op["flangeS1"] + 40), min(p[0] for p in starts) - 15)   # `starts`: last panel's start points
axs[1].set_xlim(xmin, i["cavityR"] + o["wallT"] + 60)
axs[0].set_ylim(-(o["cavityR"] + o["wallT"] + 8), o["cavityR"] + o["wallT"] + 8)
axs[1].set_ylim(o["cavityZmin"] - o["bottomT"] - 8, o["cavityZmax"] + o["topT"] + 8)
fig.legend(handles=[
    Line2D([], [], marker="*", ls="", mfc="gold", mec="k", ms=12, label="primary vertex"),
    Line2D([], [], color="tab:purple", ls="--", label="vertex -> LXe entry (dashed only if no primary_path ntuple)"),
    Line2D([], [], marker="o", ls="", color="tab:purple", ms=4, label="first LXe step"),
    Line2D([], [], color="k", label="primary inside the LXe"),
    Line2D([], [], marker="x", ls="", color="tab:red", ms=5, label="last LXe step of the primary"),
    Line2D([], [], color="tab:green", label=f"Al window ({al['thickness']*1e3:.0f} um)"),
    Line2D([], [], color="tab:orange", label=f"Ti window ({ti['thickness']*1e3:.0f} um)"),
    Line2D([], [], color="red", lw=2.5, label="fluence planes (SiPM)")],
    loc="lower center", ncol=3, fontsize=8, framealpha=0.9)
fig.suptitle(f"beam particle paths into the cryostat, {len(events)} events "
             f"(secondaries not shown; LXe light blue, steel grey)", fontsize=11)
fig.tight_layout(rect=(0, 0.06, 1, 0.98))
fig.savefig("beam_path.png", dpi=140)
print("wrote beam_path.png;", len(events), "events;",
      "geometry: real G4 solids (geometry_map_*.txt)" if use_map else "geometry: drawn from the manifest (run ./geometryProbe --map for the real solids)", ";",
      "path: primary_path ntuple (true path)" if pp is not None else
      ("path in: straight line, vertex from tree_fuse" if vtx else f"path in: straight line, vertex assumed at x = {x0_default} mm"))
