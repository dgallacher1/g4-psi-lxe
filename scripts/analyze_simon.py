#!/usr/bin/env python3
"""Beam-profile scan (SiPM scanner, run 74): 2D Gaussian fit, comparison with the scanner's own fit and with the
default G4 beam (circular, sigma 10 mm, centred on the beam window), and the profile in G4 coordinates.

  python3 analyze_simon.py [run74_profile.csv]

CSV: row 1 = scanner fit (centre_x, centre_y, sigma_x, sigma_y [mm], norm [Hz/mA], unix time);
     rows 2.. = scan points (x, y [mm], rate density [Hz/mm^2/mA], current [mA], HV [V], I [uA]).
Rates come in steps of RATE_QUANTUM (= 1 count), so the error is Poisson: RATE_QUANTUM * sqrt(n).
The scanner picture must be rotated by 180 deg to match the real alignment (ROTATE_180).

Writes beam_profile_fit.png and prints the fit, the comparison and a GPS block for the G4 beam macro.
"""
import os
import sys

import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

# ------------------------------------------------------------------------------------------------ edit these
HERE = os.path.dirname(os.path.abspath(__file__))
CSV = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "run74_profile.csv")
RATE_QUANTUM = 0.16          # Hz/mm^2/mA per count (all rates are integer multiples of this)
ROTATE_180 = True            # scanner frame -> real alignment: (x, y) -> (-x, -y)
DROP_EMPTY_ROWS = True       # scan rows with every point exactly 0 (y = 66 mm in run 74) are treated as not measured
# Default G4 beam: circular Gaussian, sigma 10 mm, centred on the beam-window centre.
SIGMA_DEFAULT = 10.0         # mm
# Beam window = the scanner's "integration region" (circle), given in the scanner's own (unrotated) frame [mm].
# ASSUMED to be the window aperture: the default G4 beam is centred on it.
WINDOW_CENTER_RAW = (-2.5, 23.5)
WINDOW_DIAMETER = 35.0
# (rotated) scanner frame -> G4 transverse plane.  ASSUMED: G4 beam along +x, scanner vertical (PCB axis) = G4 +z,
# scanner horizontal (SiPM axis) = H_SIGN * G4 y.  The horizontal sign is the unresolved mirror ambiguity: check it.
H_SIGN = +1
# ------------------------------------------------------------------------------------------------

d = np.loadtxt(CSV, delimiter=",")
auto = d[0]
x, y, r = d[1:, 0].copy(), d[1:, 1].copy(), d[1:, 2].copy()
if DROP_EMPTY_ROWS:
    empty = np.array([np.all(r[y == v] == 0) for v in y])
    if empty.any():
        print(f"dropping {int(empty.sum())} points in empty scan rows y = {sorted(set(y[empty]))} mm")
    x, y, r = x[~empty], y[~empty], r[~empty]
n_cnt = np.round(r / RATE_QUANTUM)
err = RATE_QUANTUM * np.sqrt(np.maximum(n_cnt, 1.0))
sgn = -1.0 if ROTATE_180 else 1.0
WINDOW_CENTER_SCAN = (sgn * WINDOW_CENTER_RAW[0], sgn * WINDOW_CENTER_RAW[1])
R_WIN = WINDOW_DIAMETER / 2
X, Y = sgn * x, sgn * y                               # rotated (real-alignment) scanner frame
dx, dy = np.diff(np.unique(X)).min(), np.diff(np.unique(Y)).min()
print(f"{len(r)} scan points, grid {dx:g} x {dy:g} mm, X {X.min():g}..{X.max():g}, Y {Y.min():g}..{Y.max():g} mm "
      f"(after {'180 deg rotation' if ROTATE_180 else 'no rotation'})")


# ---------------------------------------------------------------- models
def g2(xy, A, x0, y0, s1, s2, th):
    """Rotated 2D Gaussian; s1 along the axis at angle th (rad, from +X towards +Y), s2 perpendicular."""
    X_, Y_ = xy
    c, s = np.cos(th), np.sin(th)
    u = c * (X_ - x0) + s * (Y_ - y0)
    v = -s * (X_ - x0) + c * (Y_ - y0)
    return A * np.exp(-0.5 * ((u / s1) ** 2 + (v / s2) ** 2))


def g2_axis(xy, A, x0, y0, s1, s2):
    return g2(xy, A, x0, y0, s1, s2, 0.0)


def g1(t, A, m, s):
    return A * np.exp(-0.5 * ((t - m) / s) ** 2)


def chi2(model):
    return np.sum(((r - model) / err) ** 2)


# ---------------------------------------------------------------- fits
p0 = [r.max(), np.average(X, weights=r), np.average(Y, weights=r), 25.0, 25.0, 0.0]
p_rot, c_rot = curve_fit(g2, (X, Y), r, p0=p0, sigma=err, absolute_sigma=True, maxfev=20000)
if abs(p_rot[3]) < abs(p_rot[4]):                      # convention: s1 = major axis
    p_rot[3], p_rot[4], p_rot[5] = p_rot[4], p_rot[3], p_rot[5] + np.pi / 2
    p_rot, c_rot = curve_fit(g2, (X, Y), r, p0=p_rot, sigma=err, absolute_sigma=True, maxfev=20000)
p_rot[3:5] = np.abs(p_rot[3:5])
p_rot[5] = (p_rot[5] + np.pi / 2) % np.pi - np.pi / 2   # angle in (-90, 90] deg
e_rot = np.sqrt(np.diag(c_rot))
chi_rot = chi2(g2((X, Y), *p_rot))

p_ax, c_ax = curve_fit(g2_axis, (X, Y), r, p0=p0[:5], sigma=err, absolute_sigma=True, maxfev=20000)
p_ax[3:5] = np.abs(p_ax[3:5])
e_ax = np.sqrt(np.diag(c_ax))
chi_ax = chi2(g2_axis((X, Y), *p_ax))

# 1D fits to the projections (rate summed over the other axis x grid step: Hz/mm/mA), as on the scanner plot
proj = {}
for name, A_, B_, step in (("X", X, Y, dy), ("Y", Y, X, dx)):
    u = np.unique(A_)
    pr = np.array([r[A_ == v].sum() * step for v in u])
    pe = np.array([np.sqrt(np.sum(err[A_ == v] ** 2)) * step for v in u])
    pp, _ = curve_fit(g1, u, pr, p0=[pr.max(), np.average(u, weights=pr), 25.0], sigma=pe, maxfev=20000)
    pp[2] = abs(pp[2])
    proj[name] = (u, pr, pe, pp)

npts = len(r)
norm_rot = p_rot[0] * 2 * np.pi * p_rot[3] * p_rot[4]
norm_ax = p_ax[0] * 2 * np.pi * p_ax[3] * p_ax[4]

# the scanner's fit, expressed in the rotated frame
a_x, a_y = sgn * auto[0], sgn * auto[1]
a_sx, a_sy, a_norm = auto[2], auto[3], auto[4]

# ---------------------------------------------------------------- print
print(f"\n== fits (frame: {'rotated 180 deg = real alignment' if ROTATE_180 else 'scanner'}), centres in mm ==")
print(f"{'':26}{'centre X':>10}{'centre Y':>10}{'sigma 1':>9}{'sigma 2':>9}{'angle':>8}{'norm Hz/mA':>12}{'chi2/ndf':>13}")
print(f"{'2D rotated (this)':26}{p_rot[1]:10.2f}{p_rot[2]:10.2f}{p_rot[3]:9.2f}{p_rot[4]:9.2f}"
      f"{np.degrees(p_rot[5]):7.1f}d{norm_rot:12.0f}{chi_rot:8.0f}/{npts - 6}")
print(f"{'   +-':26}{e_rot[1]:10.2f}{e_rot[2]:10.2f}{e_rot[3]:9.2f}{e_rot[4]:9.2f}{np.degrees(e_rot[5]):7.1f}d")
print(f"{'2D axis-aligned (this)':26}{p_ax[1]:10.2f}{p_ax[2]:10.2f}{p_ax[3]:9.2f}{p_ax[4]:9.2f}{'0':>7}d"
      f"{norm_ax:12.0f}{chi_ax:8.0f}/{npts - 5}")
print(f"{'1D projections (this)':26}{proj['X'][3][1]:10.2f}{proj['Y'][3][1]:10.2f}{proj['X'][3][2]:9.2f}"
      f"{proj['Y'][3][2]:9.2f}{'0':>7}d")
print(f"{'scanner auto fit (CSV)':26}{a_x:10.2f}{a_y:10.2f}{a_sx:9.2f}{a_sy:9.2f}{'0':>7}d{a_norm:12.0f}")
print(f"{'default G4 beam':26}{WINDOW_CENTER_SCAN[0]:10.2f}{WINDOW_CENTER_SCAN[1]:10.2f}{SIGMA_DEFAULT:9.2f}"
      f"{SIGMA_DEFAULT:9.2f}{'0':>7}d")
print("sigma 1/2: rotated fit = major/minor axis (angle from +X towards +Y); otherwise X/Y.")

# coverage: fraction of the fitted beam inside the scanned area
xx, yy = np.meshgrid(np.linspace(-200, 200, 801), np.linspace(-200, 200, 801))
G = g2((xx, yy), *p_rot)
inside = ((xx >= X.min() - dx / 2) & (xx <= X.max() + dx / 2) & (yy >= Y.min() - dy / 2) & (yy <= Y.max() + dy / 2))
print(f"\nfraction of the fitted beam inside the scanned area: {G[inside].sum() / G.sum():.2f}  "
      f"(sigmas larger than ~half the scan width are an extrapolation)")
off = np.array(p_rot[1:3]) - np.array(WINDOW_CENTER_SCAN)
print(f"window (integration region): D = {WINDOW_DIAMETER:g} mm at {WINDOW_CENTER_SCAN} (rotated frame)")
print(f"fitted centre relative to the window centre: ({off[0]:+.1f}, {off[1]:+.1f}) mm")
in_win = np.hypot(xx - WINDOW_CENTER_SCAN[0], yy - WINDOW_CENTER_SCAN[1]) <= R_WIN
on_grid = np.hypot(X - WINDOW_CENTER_SCAN[0], Y - WINDOW_CENTER_SCAN[1]) <= R_WIN
cell = (xx[0, 1] - xx[0, 0]) * (yy[1, 0] - yy[0, 0])
print(f"rate through the window [Hz/mA]: data (grid sum) {r[on_grid].sum() * dx * dy:.0f}, 2D fit {G[in_win].sum() * cell:.0f}; "
      f"fraction of the beam inside the window: 2D fit {G[in_win].sum() / G.sum():.2f}, "
      f"default G4 {1 - np.exp(-R_WIN ** 2 / (2 * SIGMA_DEFAULT ** 2)):.2f}")
print(f"width vs default: sigma_major {p_rot[3] / SIGMA_DEFAULT:.1f}x, sigma_minor {p_rot[4] / SIGMA_DEFAULT:.1f}x "
      f"the default {SIGMA_DEFAULT:g} mm")

# ---------------------------------------------------------------- G4 frame
y_g4, z_g4 = H_SIGN * off[0], off[1]
th_g4 = H_SIGN * p_rot[5]                             # mirror in the horizontal flips the angle
c, s = np.cos(th_g4), np.sin(th_g4)
print(f"\n== G4 frame (beam along +x; transverse y = {H_SIGN:+d} x scanner X, z = scanner Y; relative to window centre) ==")
print(f"centre (y, z) = ({y_g4:+.2f}, {z_g4:+.2f}) mm, sigma_major {p_rot[3]:.2f} mm at {np.degrees(th_g4):+.1f} deg "
      f"from +y towards +z, sigma_minor {p_rot[4]:.2f} mm")
print("GPS block (keep your macro's x of the source plane and /gps/direction; centre below is window centre + offset):")
print(f"""  /gps/pos/type Beam
  /gps/pos/centre <x_window> {y_g4:.2f} {z_g4:.2f} mm      # y, z relative to the window centre: add the window y, z
  /gps/pos/rot1 0 {c:.4f} {s:.4f}
  /gps/pos/rot2 0 {-s:.4f} {c:.4f}
  /gps/pos/sigma_x {p_rot[3]:.2f} mm
  /gps/pos/sigma_y {p_rot[4]:.2f} mm""")

# ---------------------------------------------------------------- plots
fig, axs = plt.subplots(2, 2, figsize=(13, 12))
ax = axs[0, 0]
sc = ax.scatter(X, Y, c=r, s=55, marker="s", cmap="plasma", vmin=0)
fig.colorbar(sc, ax=ax, label="rate [Hz/mm$^2$/mA]")
for k, ls in ((1, "-"), (2, "--")):
    ax.add_patch(Ellipse(p_rot[1:3], 2 * k * p_rot[3], 2 * k * p_rot[4], angle=np.degrees(p_rot[5]), fill=False,
                         ec="lime", lw=2, ls=ls, label=f"2D fit {k}$\\sigma$"))
    ax.add_patch(Ellipse((a_x, a_y), 2 * k * a_sx, 2 * k * a_sy, fill=False, ec="cyan", lw=1.5, ls=ls,
                         label=f"scanner fit {k}$\\sigma$"))
    if k == 1:
        ax.add_patch(Ellipse(WINDOW_CENTER_SCAN, WINDOW_DIAMETER, WINDOW_DIAMETER, fill=False, ec="yellow", lw=2,
                             ls="-.", label=f"window D {WINDOW_DIAMETER:g} mm"))
    ax.add_patch(Ellipse(WINDOW_CENTER_SCAN, 2 * k * SIGMA_DEFAULT, 2 * k * SIGMA_DEFAULT, fill=False, ec="w", lw=1.5,
                         ls=ls, label=f"default G4 {k}$\\sigma$"))
ax.plot(*p_rot[1:3], "+", color="lime", ms=14, mew=2)
ax.plot(a_x, a_y, "x", color="cyan", ms=10, mew=2)
ax.plot(*WINDOW_CENTER_SCAN, "+", color="w", ms=14, mew=2)
ax.add_patch(Rectangle((X.min() - dx / 2, Y.min() - dy / 2), np.ptp(X) + dx, np.ptp(Y) + dy, fill=False, ec="0.6", ls=":"))
ax.set_facecolor("0.15")
ax.set_xlim(-110, 110); ax.set_ylim(-110, 110); ax.set_aspect("equal")
ax.set_xlabel("X [mm] (SiPM axis" + (", rotated 180 deg)" if ROTATE_180 else ")"))
ax.set_ylabel("Y [mm] (PCB axis" + (", rotated 180 deg)" if ROTATE_180 else ")"))
ax.set_title("data, fits and default beam (dotted box = scanned area)")
ax.legend(fontsize=7, loc="lower left", ncol=3)

ax = axs[0, 1]
pull = (r - g2((X, Y), *p_rot)) / err
sc = ax.scatter(X, Y, c=pull, s=90, marker="s", cmap="RdBu_r", vmin=-4, vmax=4)
fig.colorbar(sc, ax=ax, label="(data - 2D fit) / error")
ax.set_aspect("equal")
ax.set_xlabel("X [mm]"); ax.set_ylabel("Y [mm]")
ax.set_title(f"pulls, rotated 2D fit: chi2/ndf = {chi_rot:.0f}/{npts - 6}")


def model_proj(f, pars, axis):
    """Projection of a 2D model summed over the scan grid (same truncation as the data)."""
    A_, B_, step = (X, Y, dy) if axis == "X" else (Y, X, dx)
    u = np.unique(A_)
    Bv = np.unique(B_)
    out = []
    for v in u:
        pts = (np.full_like(Bv, v), Bv) if axis == "X" else (Bv, np.full_like(Bv, v))
        out.append(f(pts, *pars).sum() * step)
    return u, np.array(out)


A_def = norm_rot / (2 * np.pi * SIGMA_DEFAULT ** 2)    # default shape, scaled to the fitted total rate
for ax, axis in ((axs[1, 0], "X"), (axs[1, 1], "Y")):
    u, pr, pe, pp = proj[axis]
    ax.errorbar(u, pr, yerr=pe, fmt="ko", ms=4, label="data (scan grid sum)")
    uu = np.linspace(-110, 110, 400)
    ax.plot(*model_proj(g2, p_rot, axis), "-", color="green", lw=2, label="2D rotated fit (on scan grid)")
    ax.plot(uu, g1(uu, *pp), "--", color="tab:purple", lw=1.5, label=f"1D fit: $\\mu$ {pp[1]:.1f}, $\\sigma$ {pp[2]:.1f}")
    a_m, a_s, a_o = (a_x, a_sx, a_sy) if axis == "X" else (a_y, a_sy, a_sx)
    ax.plot(uu, a_norm / (np.sqrt(2 * np.pi) * a_s) * np.exp(-0.5 * ((uu - a_m) / a_s) ** 2), ":", color="tab:cyan",
            lw=2, label=f"scanner fit (full projection): $\\mu$ {a_m:.1f}, $\\sigma$ {a_s:.1f}")
    ax.plot(*model_proj(g2, [A_def, *WINDOW_CENTER_SCAN, SIGMA_DEFAULT, SIGMA_DEFAULT, 0.0], axis), "-", color="0.5",
            lw=1.5, label=f"default G4 ($\\sigma$ {SIGMA_DEFAULT:g}), same total rate")
    ax.set_xlim(-110, 110)
    ax.set_xlabel(f"{axis} [mm]")
    ax.set_ylabel("rate [Hz/mm/mA]")
    ax.set_title(f"projection on {axis}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
fig.suptitle(f"Beam scan {os.path.basename(CSV)}: 2D Gaussian {p_rot[3]:.1f} x {p_rot[4]:.1f} mm at "
             f"{np.degrees(p_rot[5]):+.0f} deg, centre ({p_rot[1]:.1f}, {p_rot[2]:.1f}) mm")
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(CSV)), "beam_profile_fit.png")
fig.savefig(out, dpi=130)
print(f"\nwrote {out}")