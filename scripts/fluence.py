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

has_origin = "creator" in fl                              # files from before the origin columns lack it
if has_origin:
    fl["creator"] = np.asarray(fl["creator"], dtype=str)

running = {}                                              # plane name -> (u, v, weight, plane) of the weighted crossings
fig, axes = plt.subplots(len(geo["planes"]), 4, figsize=(17, 4 * len(geo["planes"])), squeeze=False)
for row, pl in enumerate(geo["planes"]):
    A = pl["area"]                                        # mm^2 (disc or rectangle, from the manifest)
    t = pl["thickness"]
    m = fl["planeID"] == pl["id"]
    ent = m & (fl["entering"] == 1)
    cos = np.abs(fl["cosTheta"][ent])
    ok = cos > 1e-3                                       # guard grazing-angle divergence
    running[pl["name"]] = (fl["u_mm"][ent][ok], fl["v_mm"][ent][ok], 1.0 / cos[ok], pl)

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
        # beam particle (parentID 0) vs shower products, and which process made the latter
        def phi(sel):                                     # crossing fluence of a subset, /mm^2/primary
            c = np.abs(fl["cosTheta"][sel]); return np.sum(1.0 / c[c > 1e-3]) / A / n_primaries
        prim = ent & (fl["parentID"] == 0)
        print(f"        primary {prim.sum():6d} crossings, fluence {phi(prim):.3e} | "
              f"secondary {(ent & ~prim).sum():6d}, fluence {phi(ent & ~prim):.3e}")
        if has_origin:
            for name in np.unique(fl["creator"][ent & ~prim]):
                sel = ent & ~prim & (fl["creator"] == name)
                print(f"          made by {name:>10}: {sel.sum():6d} crossings, fluence {phi(sel):.3e}")

    ax = axes[row]
    ax[0].hist(fl["kinE_keV"][ent], bins=100)
    ax[0].set(xlabel="entry kinetic energy [keV]", ylabel="crossings", title=f"{pl['name']} spectrum")
    ax[1].hist2d(fl["u_mm"][ent], fl["v_mm"][ent], bins=40,
                 range=[[-pl["halfU"], pl["halfU"]], [-pl["halfV"], pl["halfV"]]])
    ax[1].set(xlabel="u [mm]", ylabel="v [mm]", title=f"{pl['name']} hit map", aspect="equal")
    ax[2].hist(cos, bins=50, range=(0, 1))
    ax[2].set(xlabel="|cos(theta)| to plane normal", ylabel="crossings", title=f"{pl['name']} angles")
    if has_origin and ent.sum():                          # where were the crossing particles created?
        sec = ent & (fl["parentID"] > 0)
        ax[3].hist(fl["vtxX_mm"][sec], bins=80, log=True)
        ax[3].set(xlabel="vertex x of secondaries [mm]  (beam along +x)", ylabel="crossings",
                  title=f"{pl['name']} secondaries: creation point")
    else:
        ax[3].axis("off")

fig.tight_layout()
fig.savefig("fluence.png", dpi=130)
print("\nwrote fluence.png")

# --- running fluence vs the major dimension of a centred region -------------------------------------
# Phi(D) = (weighted crossings inside the region) / (region area) / N_primaries, where the region is the plane
# shape scaled about its centre (a disc of diameter D, or a rectangle whose long side is D, same aspect ratio).
# D -> full size gives the plane mean printed above; small D approaches the local fluence at the plane centre
# ("peak" only if the beam lands there: SiPM1 yes, SiPM2 is off the beam axis). The band is the Poisson error
# of the crossing count in the region; the curve starts at MIN_N crossings.
MIN_N, PEAK_N = 10, 100
fig2, ax2 = plt.subplots(figsize=(8.5, 5.5))
colours = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
print("\nrunning fluence (centred regions), /mm^2/primary:")
for (name, (u, v, w, pl)), col in zip(running.items(), colours):
    hu, hv = pl["halfU"], pl["halfV"]
    circ = bool(pl.get("circular", False))
    s_i = np.hypot(u, v) / hu if circ else np.maximum(np.abs(u) / hu, np.abs(v) / hv)
    order = np.argsort(s_i)
    s_sorted, cum_w = s_i[order], np.cumsum(w[order])
    d_full = 2 * hu if circ else 2 * max(hu, hv)           # diameter, or the long side
    D = np.linspace(0.02 * d_full, d_full, 400)
    s = D / d_full
    n_in = np.searchsorted(s_sorted, s, side="right")      # crossings inside the region
    phi = np.where(n_in > 0, cum_w[np.maximum(n_in - 1, 0)], np.nan) / (pl["area"] * s**2) / n_primaries
    good = n_in >= MIN_N
    if not good.any():
        print(f"  {name}: fewer than {MIN_N} crossings, skipped")
        continue
    rel = 1.0 / np.sqrt(np.maximum(n_in, 1))
    label = f"{name} ({'disc, diameter' if circ else 'rectangle, length'} {d_full:.0f} mm)"
    ax2.plot(D[good], phi[good], color=col, lw=2, label=label)
    ax2.fill_between(D[good], (phi * (1 - rel))[good], (phi * (1 + rel))[good], color=col, alpha=0.2, lw=0)
    mean = phi[-1]
    ax2.axhline(mean, color=col, ls="--", lw=0.8)
    ax2.annotate(f"mean {mean:.2e}", (D[-1], mean), textcoords="offset points", xytext=(-4, 4),
                 ha="right", fontsize=8, color=col)
    pk = np.argmax(n_in >= PEAK_N) if (n_in >= PEAK_N).any() else None
    if pk is not None:
        ax2.plot(D[pk], phi[pk], "o", color=col)
        print(f"  {name:>6}: mean {mean:.3e} | centre region D = {D[pk]:5.1f} mm ({n_in[pk]} crossings): {phi[pk]:.3e} "
              f"| centre/mean {phi[pk] / mean:.1f}")
    else:
        print(f"  {name:>6}: mean {mean:.3e} | fewer than {PEAK_N} crossings in any centred region (no stable centre value)")
ax2.set(xlabel="major dimension of the centred region: diameter (disc) or length (rectangle) [mm]",
        ylabel="mean crossing fluence in the region [/mm$^2$/primary]", yscale="log",
        title="running fluence: centred regions growing to the full plane")
ax2.legend(loc="upper right")
ax2.grid(True, which="both", alpha=0.25)
fig2.tight_layout()
fig2.savefig("fluence_running.png", dpi=130)
print("wrote fluence_running.png")
