#!/usr/bin/env python3
"""Impact of the production cut: compare the runs written by macros/cut_scan.mac.

usage: python cut_scan.py [geometry_manifest.json] [file1.root file2.root ...]
       (default files: xenon_sim_cut*um.root in the current directory)
deps : uproot numpy matplotlib      (writes cut_scan.png)

For every cut value it prints, with event-level (not per-row) errors:
  * <edep> per event in the LXe                      -> should be flat: sub-cut energy is deposited locally
  * steps per event, and the share of edep carried by steps whose particle has
    kinE < 100 keV                                   -> how much energy moved from "lumped into a parent step"
                                                        to "own low-energy track" (what a NEST yield sees)
  * crossing fluence per plane, total and secondary  -> physics seen by the SiPMs
  * everything relative to the finest cut, in sigmas
Run your NEST photon count on the same files for the number that actually matters.
"""
import glob
import json
import re
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

args = sys.argv[1:]
manifest = "geometry_manifest.json"
if args and args[0].endswith(".json"):
    manifest, args = args[0], args[1:]
geo = json.load(open(manifest))
files = args or glob.glob("xenon_sim_cut*um.root")
cutof = lambda f: int(re.search(r"cut(\d+)um", f).group(1))
files = sorted(files, key=cutof, reverse=True)                 # coarse -> fine
if len(files) < 2:
    sys.exit("need at least two xenon_sim_cut<N>um.root files (run macros/cut_scan.mac)")


def per_event_mean(values, n_ev):
    """mean and standard error of a per-event quantity (a length-n_ev array)."""
    return values.mean(), values.std(ddof=1) / np.sqrt(n_ev)


rows = []
for f in files:
    with uproot.open(f) as r:
        ev = r["events"].arrays(library="np")
        st = r["steps"].arrays(["eventID", "pdg", "edep_keV", "kinE_keV"], library="np")
        fl = r["fluence"].arrays(["eventID", "planeID", "parentID", "entering", "cosTheta"], library="np")
    n = len(ev["eventID"])
    d = {"cut": cutof(f), "n": n}
    d["edep"] = per_event_mean(ev["edep_keV"], n)
    d["steps"] = per_event_mean(ev["nSteps"].astype(float), n)
    soft = (np.abs(st["pdg"]) == 11) & (st["kinE_keV"] < 100.0)
    d["soft_share"] = per_event_mean(
        np.bincount(st["eventID"][soft], st["edep_keV"][soft], minlength=n)[:n] /
        np.maximum(ev["edep_keV"], 1e-12), n)
    for p in geo["planes"]:
        ent = (fl["planeID"] == p["id"]) & (fl["entering"] == 1) & (np.abs(fl["cosTheta"]) > 1e-3)
        w = 1.0 / np.abs(fl["cosTheta"]) / p["area"]
        for label, sel in (("all", ent), ("sec", ent & (fl["parentID"] > 0))):
            d[f"{p['name']}_{label}"] = per_event_mean(
                np.bincount(fl["eventID"][sel], w[sel], minlength=n)[:n], n)
    rows.append(d)

ref = rows[-1]                                                   # finest cut = reference
keys = ["edep", "steps", "soft_share"] + [f"{p['name']}_{l}" for p in geo["planes"] for l in ("all", "sec")]
print(f"reference: cut {ref['cut']} um ({ref['n']} events)\n")
print(f"{'quantity':>12} | " + " | ".join(f"{r['cut']:>5d} um" + " " * 8 for r in rows))
for k in keys:
    cells = []
    for r in rows:
        m, e = r[k]
        rm, re_ = ref[k]
        if r is ref or rm == 0:
            cells.append(f"{m:9.4g} (ref)  ")
        else:
            sig = (m - rm) / np.hypot(e, re_) if np.hypot(e, re_) > 0 else 0.0
            cells.append(f"{m:9.4g} {100 * (m / rm - 1):+6.2f}% {sig:+4.1f}s")
    print(f"{k:>12} | " + " | ".join(cells))
print("\n(+x.xx% = relative to the finest cut; ys = that difference in standard errors, runs treated as independent)")

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
cuts = np.array([r["cut"] for r in rows])
for k, a, title in (("edep", ax[0], "LXe deposit per event"),
                    ("soft_share", ax[1], "edep share in e± steps with kinE < 100 keV")):
    m = np.array([r[k][0] for r in rows]); e = np.array([r[k][1] for r in rows])
    a.errorbar(cuts, m, e, fmt="o-"); a.set(xscale="log", xlabel="production cut [um]", title=title)
    a.invert_xaxis()
for p in geo["planes"]:
    for l, ls in (("all", "-"), ("sec", "--")):
        m = np.array([r[f"{p['name']}_{l}"][0] for r in rows]); e = np.array([r[f"{p['name']}_{l}"][1] for r in rows])
        ax[2].errorbar(cuts, m / m[-1], e / m[-1], fmt="o" + ls, label=f"{p['name']} {l}")
ax[2].set(xscale="log", xlabel="production cut [um]", ylabel="fluence / finest-cut fluence", title="plane crossing fluence")
ax[2].invert_xaxis(); ax[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig("cut_scan.png", dpi=130)
print("wrote cut_scan.png")
