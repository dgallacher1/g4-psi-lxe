#!/usr/bin/env python3
"""Look at an Am-241 alpha-source run (macros/alpha_source.mac).

usage: python alpha_look.py [xenon_sim_alpha.root] [geometry_manifest.json]     (the manifest is accepted so that
                                                                                run.sh can call every script alike; unused)
deps : uproot numpy matplotlib      writes alpha_look.png
needs: `events` and `primary_path` (the run must have /output/path true).

Per event the primary is one alpha.  From primary_path (the primary's start point and the end point of every step,
labelled with the volume the step ran in) this gives
  * where the alphas start (should lie in the 24-25 um gold sleeve at the source position),
  * the line each alpha was drawn from (intensity check of the three GPS sources),
  * the energy with which it leaves the source into the LXe = the kinetic energy at the end of the last step before the
    first step in "LXe".  This is the quantity a vacuum alpha spectrometer sees (up to its own dead layer and resolution);
    alphas that never reach the LXe stopped in the needle.
and from `events` the energy deposited in the LXe, which should equal the exit energy (the alpha stops within ~0.1 mm).
Energies are printed positive (the spectrometer histogram is drawn as -E).
"""
import os
import re
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root_file = sys.argv[1] if len(sys.argv) > 1 else "xenon_sim_alpha.root"

# Keep in step with macros/alpha_source.mac and Geo::sourceX/Y/Z (GeometryParameters.hh).
_t = np.radians(-5.0)                                  # Geo::innerAssemblyRotZdeg
_mx, _my = -30.31, 0.0                                 # Geo::sourceModelX/Y (mm)
EXPECT = np.array([_mx * np.cos(_t) - _my * np.sin(_t),    # = Geo::sourceX/Y/Z, computed, not rounded
                   _mx * np.sin(_t) + _my * np.cos(_t),
                   0.0 - 15.0])                        # sourceModelZ + Geo::innerAssemblyShiftZ
TOL = 0.05                                             # mm, warn above this
R_IN, R_OUT = 0.0723, 0.0725                           # mm, fallback emitter radii if the GDML cannot be read
def _layer_from_gdml():
    """Emitter radii from the constants of gdml/alpha_source.gdml (R_OUT = sleeveR - tCover, R_IN = R_OUT - tAct)."""
    here = os.path.dirname(os.path.abspath(__file__))
    for path in ("gdml/alpha_source.gdml", os.path.join(here, "..", "gdml", "alpha_source.gdml")):
        try:
            c = {k: float(v) for k, v in re.findall(r'<constant\s+name="(\w+)"\s+value="([-+0-9.eE]+)"', open(path).read())}
            r_out = c["sleeveR"] - c["tCover"]
            return r_out - c["tAct"], r_out, c["tCover"], c["tAct"], path
        except (OSError, KeyError):
            continue
    return None
_g = _layer_from_gdml()
if _g:
    R_IN, R_OUT = _g[0], _g[1]
    print(f"emitter layer from {_g[4]}: tCover {1e3 * _g[2]:.3f} um, tAct {1e3 * _g[3]:.3f} um -> r = {1e3 * R_IN:.2f}..{1e3 * R_OUT:.2f} um")
else:
    print("WARNING: gdml/alpha_source.gdml not found; using the built-in emitter radii")
LINES = {5485.56: 0.848, 5442.80: 0.131, 5388.23: 0.0166}   # keV: relative intensity (as in the macro)

# The vacuum alpha-spectrometer fit (ROOT stat box), a Landau convolved with a Gaussian in x = -E [MeV]:
#   p0 = Landau width, p1 = most probable value, p2 = area (counts x MeV), p3 = Gaussian sigma.
VAC = dict(xi=0.1314, mp=-4.824, area=985.4, sig=0.105)
BIN = 0.02                      # MeV, histogram bin of the fit (the spectrometer histogram looks like 0.02 MeV)
FIT_LO, FIT_HI = -5.7, -2.0     # fit window in x = -E [MeV] (the vacuum fit runs to about -2)
MEAN_LO, MEAN_HI = 2.0, 5.7     # energy window [MeV] for the mean comparison

# --- LANDAU BEGIN ------------------------------------------------------------------------------------------------------
# Landau density without scipy: phi(l) = 1/pi * int_0^inf exp(-l t - t ln t) sin(pi t) dt, tabulated once (about 1 s) and
# interpolated.  Valid for l >= -4 (below that the density is < 1e-5 and the integral loses precision: set to 0).
LAM0 = -0.22278298              # the maximum of phi sits at l = LAM0 (ROOT's langaus.C "mpshift")
_PHI = {}


def _phi_table():
    if _PHI:
        return _PHI
    lam = np.concatenate([np.arange(-5.0, 5.0, 0.01), np.geomspace(5.0, 2000.0, 400)])
    du = 0.002
    u = np.arange(du / 2, np.sqrt(150.0), du)            # t = u^2, midpoint rule (the t ln t term is singular at t = 0)
    t = u * u
    base, w = -t * np.log(t), np.sin(np.pi * t) * 2 * u * du / np.pi
    phi = np.empty(len(lam))
    for i in range(0, len(lam), 100):
        phi[i:i + 100] = (np.exp(base[None, :] - lam[i:i + 100, None] * t[None, :]) * w[None, :]).sum(axis=1)
    phi[lam < -4.0] = 0.0
    phi = np.clip(phi, 0.0, None)
    lin = lam < 5.0
    _PHI.update(lam_lin=lam[lin], phi_lin=phi[lin], loglam=np.log(lam[~lin]), logphi=np.log(phi[~lin]),
                end=phi[-1], lam_end=lam[-1])
    return _PHI


def landau_phi(l):
    """Standard Landau density at l (maximum at LAM0, ~0.1807)."""
    p = _phi_table()
    l = np.asarray(l, float)
    out = np.zeros_like(l)
    a = (l >= -5.0) & (l < 5.0)
    out[a] = np.interp(l[a], p["lam_lin"], p["phi_lin"])
    b = (l >= 5.0) & (l <= p["lam_end"])
    out[b] = np.exp(np.interp(np.log(l[b]), p["loglam"], p["logphi"]))
    c = l > p["lam_end"]
    out[c] = p["end"] * (p["lam_end"] / l[c]) ** 2
    return out


def langaus(x, xi, mp, area, sig, n=200, ns=5.0):
    """Landau (width xi, most probable value mp) convolved with a Gaussian (sigma sig), times area: the function of
    ROOT's langaus.C, evaluated pointwise by summing over +-ns sigma."""
    x = np.asarray(x, float)
    if sig < 1e-4:
        return area * landau_phi((x - mp) / xi + LAM0) / xi
    xp = x[:, None] + sig * ns * np.linspace(-1.0, 1.0, n)[None, :]
    g = np.exp(-0.5 * ((xp - x[:, None]) / sig) ** 2)
    land = landau_phi((xp - mp) / xi + LAM0) / xi
    return area * (land * g).sum(axis=1) * (2 * ns * sig / (n - 1)) / (np.sqrt(2 * np.pi) * sig)
# --- LANDAU END --------------------------------------------------------------------------------------------------------

with uproot.open(root_file) as f:
    ev = f["events"].arrays(["eventID", "edep_keV"], library="np")
    pa = f["primary_path"].arrays(["eventID", "step", "x_mm", "y_mm", "z_mm", "kinE_keV", "volume"], library="np")

# ---- order the path rows by (event, step) -------------------------------------------------------------------------
o = np.lexsort((pa["step"], pa["eventID"]))
e = pa["eventID"][o]
step = pa["step"][o]
x, y, z, kin = (pa[k][o] for k in ("x_mm", "y_mm", "z_mm", "kinE_keV"))
vol = np.asarray(pa["volume"], dtype=str)[o]

n_ev = len(ev["eventID"])
s0 = step == 0
n_start = int(s0.sum())
x0, y0, z0, k0 = x[s0], y[s0], z[s0], kin[s0]
print(f"events {n_ev} | primary_path start points {n_start}")
if n_start != n_ev:
    print("WARNING: number of start points differs from the number of events (was /output/path true?)")

# ---- where the alphas start -----------------------------------------------------------------------------------------
mean = np.array([x0.mean(), y0.mean(), z0.mean()])
r = np.hypot(x0 - EXPECT[0], y0 - EXPECT[1])
print(f"start point mean ({mean[0]:.4f}, {mean[1]:.4f}, {mean[2]:.4f}) mm, expected ({EXPECT[0]}, {EXPECT[1]}, {EXPECT[2]})")
print(f"  radius from the needle axis: mean {1e3 * r.mean():.2f} um, min {1e3 * r.min():.2f}, max {1e3 * r.max():.2f}; "
      f"{100 * np.mean((r >= R_IN - 2e-5) & (r <= R_OUT + 2e-5)):.1f} % inside the {1e3 * R_IN:.1f}-{1e3 * R_OUT:.1f} um layer; "
      f"z {z0.min() - EXPECT[2]:+.3f} .. {z0.max() - EXPECT[2]:+.3f} mm from the needle centre")
if np.any(np.abs(mean - EXPECT) > TOL):
    print(f"WARNING: mean start point is more than {TOL} mm from the expected source position: the macro's source "
          f"position and the geometry disagree")

# ---- which line each alpha came from -------------------------------------------------------------------------------
tot = sum(LINES.values())
for el, w in LINES.items():
    print(f"  line {el:8.2f} keV: {100 * np.mean(np.abs(k0 - el) < 1.0):5.2f} % of the alphas (expected {100 * w / tot:5.2f} %)")
print(f"  started with another energy: {100 * np.mean(np.min([np.abs(k0 - el) for el in LINES], axis=0) >= 1.0):.2f} %")

# ---- energy on leaving the source: the point before the first step in the LXe ---------------------------------------
idx = np.flatnonzero(vol == "LXe")
first_ev, pos = np.unique(e[idx], return_index=True)
first_row = idx[pos]
prev = first_row - 1
ok = (first_row > 0) & (e[np.maximum(prev, 0)] == e[first_row])
exit_ev, exit_E = first_ev[ok], kin[prev[ok]]
n_exit = len(exit_E)
print(f"\nalphas reaching the LXe: {n_exit} of {n_start} ({100 * n_exit / max(n_start, 1):.1f} %); the rest stop in the needle")
if n_exit == 0:
    sys.exit("no alpha reached the LXe: check the source position, the GDML and the macro")

hb = np.arange(0, 5601, 10.0)
h = np.histogram(exit_E, hb)[0]
mode = 0.5 * (hb[np.argmax(h)] + hb[np.argmax(h) + 1])
print(f"exit energy [keV]: mean {exit_E.mean():.0f} | median {np.median(exit_E):.0f} | mode {mode:.0f} (10 keV bins) | "
      f"rms {exit_E.std():.0f} | 5 %..95 %: {np.percentile(exit_E, 5):.0f}..{np.percentile(exit_E, 95):.0f}")
for thr in (2000, 3000, 4000):
    s = exit_E > thr
    print(f"  E > {thr} keV: {100 * s.mean():5.1f} % of the escaping alphas, mean {exit_E[s].mean():.0f} keV")

# ---- energy deposited in the LXe ---------------------------------------------------------------------------------------
oi = np.argsort(ev["eventID"])
ids, dep_all = ev["eventID"][oi], ev["edep_keV"][oi]
dep = dep_all[np.searchsorted(ids, exit_ev)]
print(f"\nLXe deposit of the escaping alphas [keV]: mean {dep.mean():.0f} | mean deposit / exit energy {np.mean(dep / exit_E):.4f} "
      f"| deposit < 99 % of the exit energy in {100 * np.mean(dep < 0.99 * exit_E):.2f} % of them")
print(f"events with an LXe deposit: {int(np.sum(dep_all > 0))} (alphas reaching the LXe: {n_exit})")
nodep = dep <= 0
if nodep.any():
    print(f"  {int(nodep.sum())} alphas reach the LXe without any deposit; their exit energy [keV]: "
          f"median {np.median(exit_E[nodep]):.3g}, max {exit_E[nodep].max():.3g}")

# ---- figure --------------------------------------------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
ax[0, 0].hist(exit_E, bins=hb[::2], histtype="stepfilled", alpha=0.8)
for el in LINES:
    ax[0, 0].axvline(el, color="k", lw=0.6, ls="--")
ax[0, 0].set(xlabel="energy on leaving the source [keV]", ylabel="alphas", yscale="log",
             title=f"exit energy of {n_exit} escaping alphas (dashed: Am-241 lines)")
ax[0, 1].hist(dep, bins=120, histtype="stepfilled", alpha=0.8)
ax[0, 1].set(xlabel="energy deposited in the LXe [keV]", ylabel="alphas", yscale="log", title="LXe deposit per escaping alpha")
ax[1, 0].hist(1e3 * r, bins=60)
ax[1, 0].axvline(1e3 * R_IN, color="k", ls="--", lw=0.8)
ax[1, 0].axvline(1e3 * R_OUT, color="k", ls="--", lw=0.8)
ax[1, 0].set(xlabel="start radius from the needle axis [um]", ylabel="alphas", title="where the alphas start (dashed: layer)")
ax[1, 1].hist(z0 - EXPECT[2], bins=60)
ax[1, 1].set(xlabel="start z relative to the needle centre [mm]", ylabel="alphas", title="start position along the needle")
fig.tight_layout()
fig.savefig("alpha_look.png", dpi=130)
print("\nwrote alpha_look.png")

# ---- Landau (x) Gaussian fit of the LXe deposit, compared with the vacuum fit ---------------------------------------------
# Same function and variable as the vacuum fit: x = -E [MeV], histogram bins of BIN, chi-square fit with sqrt(N) errors on
# the non-empty bins of [FIT_LO, FIT_HI].  The sim has no detector, so its Gaussian sigma is the width the fit needs to
# describe the line structure and the thin-layer edge, not a resolution.
try:
    from scipy.optimize import curve_fit
except ImportError:
    curve_fit = None
    print("scipy not available: no Landau fit")

if curve_fit is not None:
    x_dep = -dep_all[dep_all > 0] / 1e3                       # all events with an LXe deposit, x = -E in MeV
    if len(x_dep) < 3000:
        print(f"\nNOTE: only {len(x_dep)} entries for the Landau fit; use the 100k run for a stable fit")
    edges = np.arange(-6.0, BIN / 2, BIN)
    cnt, edges = np.histogram(x_dep, edges)
    cen = 0.5 * (edges[1:] + edges[:-1])
    win = (cen >= FIT_LO) & (cen <= FIT_HI)
    sel = win & (cnt > 0)
    try:
        pk = cen[win][np.argmax(cnt[win])]
        p0 = [0.1, pk, cnt[win].max() * 0.1 / 0.18, 0.05]
        popt, pcov = curve_fit(langaus, cen[sel], cnt[sel], p0=p0, sigma=np.sqrt(cnt[sel]), absolute_sigma=True,
                               bounds=([0.005, -6.0, 0.0, 0.0], [1.0, -2.0, np.inf, 0.5]), max_nfev=2000)
        perr = np.sqrt(np.diag(pcov))
        chi2 = np.sum(((cnt[sel] - langaus(cen[sel], *popt)) / np.sqrt(cnt[sel])) ** 2)
        ndf = int(sel.sum()) - 4
        sim = dict(xi=popt[0], mp=popt[1], area=popt[2], sig=popt[3])

        def wmean(par):                                        # mean energy [MeV] of a Langaus in the energy window
            xg = np.linspace(-MEAN_HI, -MEAN_LO, 3000)
            fv = langaus(xg, par["xi"], par["mp"], par["area"], par["sig"])
            return -(xg * fv).sum() / fv.sum()

        e_dep = -x_dep
        in_w = (e_dep >= MEAN_LO) & (e_dep <= MEAN_HI)
        m_data, m_sim, m_vac = e_dep[in_w].mean(), wmean(sim), wmean(VAC)
        print(f"\n=== Landau (x) Gaussian fit of the LXe deposit (x = -E), window [{FIT_LO}, {FIT_HI}] MeV, {int(sel.sum())} bins ===")
        print(f"chi2/ndf = {chi2:.1f} / {ndf}")
        print(f"{'':28s}{'LXe deposit (sim)':>26s}{'vacuum fit':>14s}")
        print(f"{'Landau width  p0 [MeV]':28s}{popt[0]:>18.4f} +- {perr[0]:.4f}{VAC['xi']:>14.4f}")
        print(f"{'most probable x  p1 [MeV]':28s}{popt[1]:>18.4f} +- {perr[1]:.4f}{VAC['mp']:>14.4f}")
        print(f"{'Gaussian sigma  p3 [MeV]':28s}{popt[3]:>18.4f} +- {perr[3]:.4f}{VAC['sig']:>14.4f}")
        print(f"{'area p2 (counts x MeV)':28s}{popt[2]:>18.1f} +- {perr[2]:.1f}{VAC['area']:>14.1f}   (depends on N and the bin: not comparable)")
        print(f"most probable energy: sim {-popt[1]:.3f} MeV, vacuum {-VAC['mp']:.3f} MeV: the vacuum source peaks "
              f"{(-popt[1]) - (-VAC['mp']):+.3f} MeV lower ({(-VAC['mp']) / (-popt[1]):.3f} x the sim)")
        print(f"mean energy in [{MEAN_LO}, {MEAN_HI}] MeV: sim data {m_data:.3f} | sim fit {m_sim:.3f} | vacuum fit {m_vac:.3f}  "
              f"-> vacuum / sim = {m_vac / m_data:.3f} (sim data), {m_vac / m_sim:.3f} (sim fit)")

        xf = np.linspace(-6.0, -1.0, 600)
        vac_c = langaus(cen[win], VAC["xi"], VAC["mp"], VAC["area"], VAC["sig"])
        scale = cnt[win].sum() / vac_c.sum()                   # same number of counts in the fit window
        fig2, ax2 = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw=dict(height_ratios=[3, 1]), sharex=True)
        ax2[0].step(cen, cnt, where="mid", color="tab:blue", label=f"sim: LXe deposit, {len(x_dep)} alphas")
        ax2[0].plot(xf, langaus(xf, *popt), "r", lw=2, label=f"sim Landau (x) Gauss fit, chi2/ndf {chi2:.0f}/{ndf}")
        ax2[0].plot(xf, scale * langaus(xf, VAC["xi"], VAC["mp"], VAC["area"], VAC["sig"]), "k--", lw=2,
                    label="vacuum spectrometer fit (shape, scaled to the same counts in the window)")
        ax2[0].axvspan(FIT_LO, FIT_HI, color="0.9", zorder=0)
        ax2[0].set(ylabel=f"alphas per {BIN} MeV", title="LXe deposit of the Am-241 alphas against the vacuum Landau fit")
        ax2[0].legend(loc="center right", fontsize=9)
        ax2[0].text(0.98, 0.95, f"sim:    xi {popt[0]:.3f}  MPV {popt[1]:.3f}  sigma {popt[3]:.3f}\n"
                                f"vacuum: xi {VAC['xi']:.3f}  MPV {VAC['mp']:.3f}  sigma {VAC['sig']:.3f}",
                    transform=ax2[0].transAxes, ha="right", va="top", family="monospace", fontsize=9)
        fit_c = langaus(cen, *popt)
        ok = cnt > 0
        ax2[1].axhline(0, color="k", lw=0.6)
        ax2[1].plot(cen[ok & win], ((cnt - fit_c) / np.sqrt(np.maximum(cnt, 1)))[ok & win], "r.", ms=4)
        ax2[1].set(xlabel="-E [MeV]", ylabel="(data - fit) / error", xlim=(-6.0, -1.0))
        fig2.tight_layout()
        fig2.savefig("alpha_fit.png", dpi=130)
        print("wrote alpha_fit.png")
    except Exception as err:                                   # a failed fit must not hide the results above
        print(f"Landau fit failed: {err}")
