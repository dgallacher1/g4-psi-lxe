"""Cross-check steps / tree_chroma / tree_fuse of one output file.   python scripts/smoke_check.py xenon_sim_smoke.root"""
import sys, numpy as np, pandas as pd, awkward as ak, uproot


def table(tree, cols):
    """Flat ntuple -> DataFrame without library="pd" (which needs awkward-pandas)."""
    a = tree.arrays(cols, library="ak")
    return pd.DataFrame({c: ak.to_numpy(ak.flatten(a[c], axis=None)) for c in cols})


f = uproot.open(sys.argv[1] if len(sys.argv) > 1 else "xenon_sim_smoke.root")
print("trees:", [k for k in f.keys()])
ev = f["events"].arrays(library="np")
e_ev = np.zeros(ev["eventID"].max() + 1); e_ev[ev["eventID"]] = ev["edep_keV"]

ch = table(f["tree_chroma"], ["event", "hits", "edep"])
fz = table(f["tree_fuse"], ["eventid", "trackid", "parentid", "type", "parenttype", "creaproc", "edproc", "ed",
                            "xp", "x", "xp_pri", "x_pri"])
for c in ("type", "parenttype", "creaproc", "edproc"):
    fz[c] = fz[c].astype(str)

s_ch = ch.groupby("event")["edep"].sum().reindex(range(len(e_ev)), fill_value=0).values
s_fz = fz.groupby("eventid")["ed"].sum().reindex(range(len(e_ev)), fill_value=0).values
print(f"rows: chroma {len(ch)}  fuse {len(fz)}  (fuse - chroma = zero-edep gamma rows: {len(fz)-len(ch)})")
print("max |events - chroma| keV:", np.abs(e_ev - s_ch).max(), "  max |events - fuse| keV:", np.abs(e_ev - s_fz).max())
print("fuse rows with ed==0:", (fz.ed == 0).sum(), " all of them gamma:", bool((fz[fz.ed == 0].type == "gamma").all()))
print("type 'unknown':", (fz.type == "unknown").sum(), "(expect 0)")
print("hits per event start at 0 and are consecutive:",
      bool(ch.groupby("event")["hits"].apply(lambda h: sorted(h) == list(range(len(h)))).all()))
print("xp == x, xp_pri == x_pri:", bool(np.allclose(fz.xp, fz.x)), bool(np.allclose(fz.xp_pri, fz.x_pri)),
      "  primary vertex x:", fz.xp_pri.unique()[:3])
print("\ncreator process x type (edep keV):")
print(fz.groupby(["creaproc", "type"])["ed"].agg(["count", "sum"]).sort_values("sum", ascending=False).head(15).to_string())
print("\nedproc counts:\n", fz.edproc.value_counts().head(10).to_string())
print("\nparenttype of primaries (parentid==0):", fz[fz.parentid == 0].parenttype.unique(),
      fz[fz.parentid == 0].creaproc.unique())

if "run_info" in f:                                       # provenance written by RunMeta
    ri = table(f["run_info"], ["key", "value"])
    print("\nrun_info:")
    for k, v in zip(ri.key, ri.value):
        if not str(k).startswith("cmd_"):
            print(f"  {k:<46} {v}")
    cmds = [str(v) for k, v in zip(ri.key, ri.value) if str(k).startswith("cmd_")]
    print(f"  commands recorded: {len(cmds)}; last 8:")
    for c in cmds[-8:]:
        print("    ", c)
else:
    print("\nno run_info ntuple (file from an older build)")
