#!/usr/bin/env python3
"""Overlay a few G4 beam events (tree_chroma energy deposits) on the chroma geometry.   python look_at_tracks.py

Edit the constants, run, open look_tracks.png.  Geometry loading is the same as look_at_geometry.py (Absolute Zero on the
Center component, components from the YAML list), so keep the two files together.

G4 -> chroma map:  p_chroma = A @ p_g4 + B.   Start with identity and set B[2] (dz) / flip signs in A until the beam enters
through the window in the right place.  G4: beam along +x from x = -175 mm, z = 0 is the beam-window centre.
"""
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import look_at_geometry as G

# ------------------------------------------------------------------------------------------------ edit these
G4_ROOT = "xenon_sim_focused.root"        # G4 output with tree_chroma
N_EVENTS = 5                              # events to draw (the first ones with deposits)
A = np.eye(3)                             # G4 -> chroma:  p_chroma = A @ p_g4 + B     (e.g. np.diag([-1, -1, 1]))
B = np.array([0.0, 0.0, 0.0])             #                B[2] = the z offset
OUT = "look_tracks.png"
# ------------------------------------------------------------------------------------------------


def load_events():
    f = uproot.open(G4_ROOT)
    name = next((k for k in f.keys() if k.split(";")[0] == "tree_chroma"), None)
    if name is None:
        raise SystemExit(f"{G4_ROOT} has no tree_chroma (keys: {f.keys()})")
    t = f[name]
    n = min(t.num_entries, 5000)
    a = t.arrays(["event", "xpos", "ypos", "zpos", "edep"], entry_stop=n, library="np")
    ev = []
    for i in range(n):
        x, y, z, e = (np.atleast_1d(np.asarray(a[k][i], float)) for k in ("xpos", "ypos", "zpos", "edep"))
        if len(x) and e.sum() > 0:
            ev.append((int(np.atleast_1d(a["event"][i])[0]), np.stack([x, y, z], axis=1), e))
        if len(ev) == N_EVENTS:
            break
    if not ev:
        raise SystemExit("no events with deposits in the first 5000 entries")
    return ev


def main():
    c = G.absolute_zero()
    print(f"Absolute Zero c = {np.array2string(c, precision=5)}")
    comps = {n: [(f, G.read_stl(f) + c) for f in G.find_files(n)] for n in G.COMPONENTS}
    cols = [plt.get_cmap("tab20")(i % 20) for i in range(len(comps))]
    rng = np.random.default_rng(1)
    events = load_events()
    ecol = plt.get_cmap("Set1")

    fig, axs = plt.subplots(2, 2, figsize=(16, 14))
    for ax, view in zip((axs[0, 0], axs[0, 1], axs[1, 0]), ("x-z", "x-y", "y-z")):
        ia, ib = G.VIEWS[view]
        for (name, items), col in zip(comps.items(), cols):
            if items:
                v = np.concatenate([t.reshape(-1, 3) for _, t in items])
                v = v[rng.choice(len(v), min(len(v), G.MAX_PTS), replace=False)]
                ax.scatter(v[:, ia], v[:, ib], s=0.4, color=col, alpha=0.35, linewidths=0)
        for k, (eid, p, e) in enumerate(events):
            q = p @ A.T + B
            ax.plot(q[:, ia], q[:, ib], "-", lw=0.4, color=ecol(k % 9), alpha=0.7)
            ax.scatter(q[:, ia], q[:, ib], s=4 + 30 * e / e.max(), color=ecol(k % 9), alpha=0.8, linewidths=0, label=f"event {eid}" if view == "x-z" else None)
        G.source_cylinder(ax, view, G.SOURCE_CFG, color="k", lw=1.0)
        G.decorate(ax, view)
        ax.set_title(f"{view}: chroma frame; coloured = G4 deposits (marker size ~ edep), black = configured source", fontsize=9)
    axs[0, 0].legend(fontsize=7, markerscale=1.5)
    ax = axs[1, 1]
    ax.axis("off")
    lines = [f"p_chroma = A p_g4 + B,  B = {B}", f"A =\n{A}", ""]
    for eid, p, e in events:
        q = p @ A.T + B
        cen = (q * e[:, None]).sum(axis=0) / e.sum()
        lines.append(f"event {eid}: {len(p)} deposits, sum {e.sum():.0f} keV, x {q[:, 0].min():.1f}..{q[:, 0].max():.1f}, "
                     f"z {q[:, 2].min():.1f}..{q[:, 2].max():.1f}, edep centroid ({cen[0]:.1f}, {cen[1]:.1f}, {cen[2]:.1f}) mm")
    ax.text(0.0, 1.0, "\n".join(lines), fontsize=8, va="top", family="monospace", transform=ax.transAxes)
    fig.tight_layout()
    fig.savefig(OUT, dpi=110)
    print("\n".join(lines))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
