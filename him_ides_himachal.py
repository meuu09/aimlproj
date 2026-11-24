"""
Temporal evolution of seismicity in the Himachal-Kashmir Himalaya (USGS ComCat)
Adapted from kemichai/Him_ides (Michailos et al. 2021, Frontiers in Earth Science).

What is kept from the original code (python/temporal_plots.py, Figures_5-9.py):
  * magnitude-vs-time + cumulative count plot with a locator map   (orig. Fig. 5)
  * inter-event time plot (log scale) with a locator map           (orig. Fig. 6)
  * magnitude-frequency distribution + b-value from a straight-line
    fit to log10(cumulative number), using max-curvature for Mc    (orig. Fig. 7)
  * migration plot: distance along the arc vs time                 (orig. Figs. 8-9, simplified)
  * 3-D hypocentre plot                                            (orig. Fig. 10, simplified)

What is changed:
  * input is a USGS CSV instead of a QuakeML file (no ObsPy Catalog needed)
  * dates, map extents, region box and magnitude cutoff are parameters below
  * Mc by maximum curvature is computed two ways, so eqcorrscan does not have to
    be installed: (1) the textbook Wiemer & Wyss (2000) version, the most
    frequent 0.1-magnitude bin, and (2) the exact eqcorrscan 0.4.2 calc_max_curv
    code the authors used (peak of the second derivative on 0.5 bins), with the
    paper's +0.4 correction (Woessner & Wiemer 2005)
  * only events whose depth USGS actually solved for are used (origin
    "depth-type" flag in data/usgs_depth_type.csv, from fetch_depth_type.py).
    Fixed ("operator assigned") and unflagged depths are dropped entirely.
  * an event counts as intermediate only if it stays below 40 km after
    subtracting its depth uncertainty: depth - depthError >= 40 km
  * migration uses a profile fitted to the data (PCA) instead of hand-picked lines
  * extra figure: rate / Mc / b-value per decade (the temporal-evolution summary)

Run from anywhere:  .venv/bin/python him_ides_himachal.py                   (results/excluding_fixed/)
                    add --catalog combined for USGS + ISC   (results/combined_usgs_isc/<case>/)
                    .venv/bin/python him_ides_himachal.py --include-fixed   (results/including_fixed/)
                    --catalog combined --drop-duplicates  (results/combined_usgs_isc_nodup/<case>/)
Windows (PowerShell): .venv\\Scripts\\python him_ides_himachal.py [same options]
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from seis_utils import (BOX, CATALOG_END, DEEP_MAX, DEEP_MIN, LONG_CLASS, SOLVED_DEPTH, b_lsq, b_mle, case_label,
                        case_name, cv, interevent_days, load_catalog, max_curv, max_curv_eqcorrscan,
                        parse_args, results_base)

# ----------------------------- PARAMETERS ---------------------------------- #
# Options (see seis_utils.parse_args):
#   (default)          excluding_fixed: only solved depths; fixed/unflagged events dropped
#   --include-fixed    including_fixed: every event, depths as reported
#   --catalog combined USGS + ISC Bulletin (data/combined_usgs_isc_catalog.csv)
#   --drop-duplicates  combined only: drop rows flagged as possible duplicates
# Both depth cases use the same intermediate rule (depth - depthError >= 40 km).
ARGS = parse_args(__doc__.splitlines()[1])
CATALOG, INCLUDE_FIXED = ARGS.catalog, ARGS.include_fixed
CASE = case_name(ARGS)
USED = "all depths" if INCLUDE_FIXED else "solved depths"   # wording for labels
CASE_LABEL = case_label(ARGS)
OUT = os.path.join(results_base(ARGS), CASE)
START_YEAR = 1990          # analysis period
MC_CUT = 4.4               # magnitude cutoff used for rate / inter-event / b-value
MAG_WORD = "mostly mb" if CATALOG == "usgs" else "mixed types: mb, ML, MD, Mw"
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    """Save a figure into this case's folder, stamped with the case name."""
    fig.text(0.005, 0.995, CASE_LABEL, ha="left", va="top", fontsize=8, color="0.4")
    fig.savefig(os.path.join(OUT, name), dpi=200)

C_DEEP, C_SHAL, C_FIX = "#d62728", "#1f77b4", "#9a9a9a"
matplotlib.rc("font", size=13)

# ----------------------------- LOAD / PREPARE ------------------------------ #
raw, load_info = load_catalog(CATALOG, INCLUDE_FIXED, START_YEAR, ARGS.drop_duplicates)
n_before = load_info["after_duplicates"]
print(f"{CASE_LABEL} | {START_YEAR}+: {n_before} events, {len(raw)} used, "
      f"{n_before - len(raw)} fixed/unflagged dropped")
# depth_class: short names from seis_utils, written out with their long labels
raw["depth_class"] = raw.depth_class.map(LONG_CLASS)
MAG_MIN = raw.mag.min()                     # lowest magnitude actually present (labels)
cat = raw[raw.mag >= MC_CUT].sort_values("time").reset_index(drop=True)
deep = cat[cat.depth_class.str.startswith("intermediate")]
shal = cat[cat.depth_class.str.startswith("shallow")]
unc = cat[cat.depth_class.str.startswith("uncertain")]
vdeep = cat[cat.depth_class.str.startswith("deep")]          # > 150 km
print(f"Catalog {START_YEAR}+, M>={MC_CUT}: N={len(cat)}  intermediate={len(deep)}  "
      f"shallow={len(shal)}  uncertain={len(unc)}  deeper than {DEEP_MAX:.0f} km={len(vdeep)}")
print(cat.depth_class.value_counts().to_string())
print("Magnitude types, M>=%.1f:" % MC_CUT, cat.magType.value_counts().to_dict())


def locator(ax, df_groups):
    for df, col, mk in df_groups:
        ax.scatter(df.longitude, df.latitude, s=14, c=col, marker=mk, alpha=0.8, lw=0.6)
    ax.set_xlim(BOX["lon"][0] - 0.3, BOX["lon"][1] + 0.3)
    ax.set_ylim(BOX["lat"][0] - 0.3, BOX["lat"][1] + 0.3)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.yaxis.set_label_position("right"); ax.yaxis.tick_right()
    ax.tick_params(top=True, left=True)


# max_curv, max_curv_eqcorrscan, b_lsq, b_mle: see seis_utils.py

# ----------------- FIG A: magnitude vs time + cumulative ------------------- #
fig = plt.figure(figsize=(13, 6.5))
ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=3, rowspan=2)
groups = [(unc, C_FIX, "+", f"Uncertain: error bar crosses 40 km (N={len(unc)})"),
          (shal, C_SHAL, "o", f"Shallow: depth + error < 40 km (N={len(shal)})"),
          (deep, C_DEEP, "X", f"Intermediate: depth − error ≥ 40 km (N={len(deep)})")]
if len(vdeep):
    groups.append((vdeep, "k", "D", f"Deeper than {DEEP_MAX:.0f} km (N={len(vdeep)})"))
for df, col, mk, lab in groups:
    ax1.vlines(df.time, 0, df.mag, color="gray", lw=0.3, zorder=1)
    ax1.scatter(df.time, df.mag, s=df.mag**2 * 2.2, c=col if mk != "o" else "none",
                edgecolors=col, marker=mk, alpha=0.8, lw=1.4, zorder=3, label=lab)
ax1.set_ylim(MC_CUT - 0.2, cat.mag.max() + 0.3)
ax1.set_xlim(pd.Timestamp(f"{START_YEAR}-01-01"), cat.time.max() + pd.Timedelta(days=200))
ax1.set_ylabel(f"Magnitude ({MAG_WORD})"); ax1.set_xlabel("Time (UTC)")
ax1.xaxis.set_major_locator(mdates.YearLocator(5)); ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax1.legend(loc="upper left", fontsize=10)
ax2 = ax1.twinx()
ax2.plot(cat.time, np.arange(1, len(cat) + 1), "k-", lw=2.2, label="All, M≥%.1f" % MC_CUT)
ax2.plot(deep.time, np.arange(1, len(deep) + 1), "k--", lw=2.2, label="Intermediate depth")
ax2.set_ylabel("Cumulative number", color="gray"); ax2.tick_params(axis="y", labelcolor="gray")
ax2.legend(loc="lower right", fontsize=10)
ax4 = plt.subplot2grid((2, 4), (0, 3))
locator(ax4, [(unc, C_FIX, "+"), (shal, C_SHAL, "o"), (deep, C_DEEP, "X")])
plt.tight_layout(); save(fig, "Fig1_magnitude_vs_time_cumulative.png"); plt.close(fig)

# --------------------- FIG B: inter-event times ---------------------------- #
def interevent(df):
    return interevent_days(df.time)

fig = plt.figure(figsize=(13, 6.5))
ax = plt.subplot2grid((2, 4), (0, 0), colspan=3, rowspan=2)
ta, da = interevent(cat); td, dd = interevent(deep)
ax.scatter(ta, da, s=45, marker="+", c=C_FIX, lw=1.2, label=f"All events ({USED}), M≥{MC_CUT} (N={len(cat)})", zorder=1)
ax.scatter(td, dd, s=70, marker="o", facecolors="none", edgecolors=C_DEEP, lw=1.8,
           label=f"Intermediate depth only (N={len(deep)})", zorder=3)
ax.set_yscale("log"); ax.set_ylabel("Inter-event time (days)"); ax.set_xlabel("Time (UTC)")
ax.set_xlim(pd.Timestamp(f"{START_YEAR}-01-01"), cat.time.max() + pd.Timedelta(days=200))
ax.xaxis.set_major_locator(mdates.YearLocator(5)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.legend(loc="lower left", fontsize=10); ax.tick_params(top=True, right=True)
ax4 = plt.subplot2grid((2, 4), (0, 3))
locator(ax4, [(cat, C_FIX, "+"), (deep, C_DEEP, "o")])
plt.tight_layout(); save(fig, "Fig2_interevent_times.png"); plt.close(fig)
print(f"Inter-event time, all: median={np.median(da):.1f} d, mean={da.mean():.1f} d, "
      f"min={da.min():.3f} d, max={da.max():.0f} d | CV={cv(da):.2f}")
print(f"Inter-event time, intermediate: median={np.median(dd):.1f} d, mean={dd.mean():.1f} d")

# ---------------------- FIG C: magnitude-frequency + b --------------------- #
mags_all = raw.mag.values                      # every event used in this case (1990+) to locate Mc
mc_curv = max_curv(mags_all)
fig, ax1 = plt.subplots(figsize=(9.5, 6.5))
bins = np.arange(np.floor(MAG_MIN * 10) / 10 - 0.05, raw.mag.max() + 0.2, 0.2)
ax1.hist(raw.mag, bins, histtype="step", fill=True, color="gray", alpha=0.25, lw=2, ec="k",
         label=f"All events ({USED}), M≥{MAG_MIN:.1f} since {START_YEAR} (N={len(raw)})")
ax1.set_ylabel("Frequency"); ax1.set_xlabel(f"Magnitude ({MAG_WORD})")
ax2 = ax1.twinx()
res = {}
for lab, df, col, ls, mk in [(f"All events ({USED})", raw, "k", "-", "+"),
                             ("Intermediate (depth − error ≥ 40 km)",
                              raw[raw.depth_class.str.startswith("intermediate")], C_DEEP, "--", "o")]:
    m = df.mag.values
    b, steps, cum, fn = b_lsq(m, MC_CUT)
    bm, bs, n = b_mle(m, MC_CUT)
    res[lab] = (b, bm, bs, n)
    allsteps = np.unique(m)
    allcum = np.array([(m >= s).sum() for s in allsteps])
    ax2.scatter(allsteps, np.log10(allcum), marker=mk, s=40, c=col if mk == "+" else "none",
                edgecolors=col, alpha=0.7, lw=1.3)
    if fn is not None:
        ax2.plot(steps, fn(steps), ls, c=col, lw=2.4,
                 label=f"{lab}: b={b:.2f} (fit), {bm:.2f}±{bs:.2f} (MLE)\nM≥{MC_CUT}, N={n}")
ax2.axvline(MC_CUT, color="k", lw=0.8, ls=":")
ax2.set_ylabel("log$_{10}$ cumulative number")
# legend below the axes: inside the plot it covered the cumulative curve
ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=10, ncol=1, frameon=False)
ax1.set_title(f"Max-curvature Mc of this catalog = {mc_curv}; fit threshold M≥{MC_CUT}", fontsize=11)
plt.tight_layout(); save(fig, "Fig3_magnitude_frequency_bvalue.png"); plt.close(fig)
print(f"Max-curvature Mc (all used events, {START_YEAR}+): Wiemer-Wyss mode = {mc_curv}; "
      f"authors' eqcorrscan method = {max_curv_eqcorrscan(mags_all)} (+0.4 = {max_curv_eqcorrscan(mags_all) + 0.4:.1f})")
for k, v in res.items():
    print(f"b-value [{k}]: LSQ={v[0]:.2f}  MLE={v[1]:.2f}±{v[2]:.2f}  N={v[3]}")

# ------------------------ FIG D: migration along arc ----------------------- #
lat0, lon0 = cat.latitude.mean(), cat.longitude.mean()
x = (cat.longitude - lon0) * 111.32 * np.cos(np.radians(lat0))
y = (cat.latitude - lat0) * 110.57
u, s, vt = np.linalg.svd(np.c_[x - x.mean(), y - y.mean()], full_matrices=False)
ax_dir = vt[0] if vt[0][1] < 0 else -vt[0]          # point toward SE so distance grows SE
cat["along_km"] = np.c_[x, y] @ ax_dir
strike_az = (np.degrees(np.arctan2(ax_dir[0], ax_dir[1])) + 360) % 360
fig, ax = plt.subplots(figsize=(12, 6))
sc = ax.scatter(cat.time, cat.along_km, s=cat.mag**2 * 2.5,
                c=cat.depth, cmap="viridis_r", vmin=0, vmax=150,
                edgecolors=np.where(cat.depth_class.str.startswith("intermediate"), C_DEEP, "k"),
                linewidths=np.where(cat.depth_class.str.startswith("intermediate"), 2.0, 0.4),
                alpha=0.85, zorder=3)
ax.scatter([], [], s=60, facecolors="none", edgecolors=C_DEEP, lw=2, label="intermediate (depth − error ≥ 40 km)")
cb = plt.colorbar(sc, ax=ax)
cb.set_label("Solved depth (km)" if not INCLUDE_FIXED else "Depth as reported, incl. fixed (km)")
ax.set_xlabel("Time (UTC)"); ax.set_ylabel(f"Distance along arc, azimuth {strike_az:.0f}° (km)")
ax.xaxis.set_major_locator(mdates.YearLocator(5)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.legend(loc="upper left", fontsize=10); ax.grid(alpha=0.2)
plt.tight_layout(); save(fig, "Fig4_migration_along_arc.png"); plt.close(fig)
print(f"Profile azimuth (PCA of epicentres): {strike_az:.0f} deg; along-arc extent "
      f"{cat.along_km.min():.0f} to {cat.along_km.max():.0f} km")

# ----------------------------- FIG E: 3-D ---------------------------------- #
meas = cat
fig = plt.figure(figsize=(10, 7.5))
ax = fig.add_subplot(111, projection="3d")
p = ax.scatter(meas.longitude, meas.latitude, -meas.depth, c=meas.depth, cmap="viridis_r", vmin=0, vmax=150,
               s=meas.mag**2 * 2.5, edgecolors="k", linewidths=0.3, alpha=0.9)
ax.set_xlabel("Longitude (°)"); ax.set_ylabel("Latitude (°)"); ax.set_zlabel("Depth (km)")
ax.set_title((f"Solved depths only, M≥{MC_CUT} (N={len(meas)}); fixed/unflagged depths not used"
              if not INCLUDE_FIXED else
              f"All depths as reported, M≥{MC_CUT} (N={len(meas)}); includes fixed/unflagged depths"),
             fontsize=11)
fig.colorbar(p, ax=ax, shrink=0.6, label="Depth (km)"); ax.view_init(elev=22, azim=-65)
plt.tight_layout(); save(fig, "Fig5_3D_hypocentres.png"); plt.close(fig)

# ------------- FIG F: temporal evolution summary (rate, Mc, b) ------------- #
# CAUTION (see results/qc/QC_REPORT.md): rates per period mix real seismicity with
# changes in detection and in how often depths were solved. Do not read trends from them.
periods = [(1990, 1999), (2000, 2009), (2010, 2019), (2020, 2026)]
t_end = CATALOG_END                     # end of the downloads, not the last event kept
rows = []
for a, b_ in periods:
    seg = raw[(raw.time.dt.year >= a) & (raw.time.dt.year <= b_)]
    # period length in years; the last period ends at the last event in the catalog
    stop = min(pd.Timestamp(f"{b_ + 1}-01-01"), t_end)
    yrs = (stop - pd.Timestamp(f"{a}-01-01")).days / 365.25
    n44 = (seg.mag >= MC_CUT).sum()
    bm, bs, n = b_mle(seg.mag.values, MC_CUT)
    rows.append(dict(period=f"{a}-{b_}", Mc_maxcurv=max_curv(seg.mag.values),
                     N_all_used=len(seg), mag_min=seg.mag.min(), N_Mc=n44,
                     share_solved_depth_Mc=(seg[seg.mag >= MC_CUT].depth_type.isin(SOLVED_DEPTH).mean()
                                            if n44 else np.nan),
                     years=round(yrs, 2), rate_per_yr=n44 / yrs, b_MLE=bm, b_err=bs))
tab = pd.DataFrame(rows); tab.to_csv(os.path.join(OUT, "per_period_summary.csv"), index=False)
print(tab.round(2).to_string(index=False))
fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
xs = np.arange(len(tab))
axs[0].bar(xs, tab.rate_per_yr, color="#555"); axs[0].set_ylabel(f"Events per year (M≥{MC_CUT})")
axs[1].plot(xs, tab.Mc_maxcurv, "o-", c="k"); axs[1].set_ylabel("Mc (max curvature)")
axs[1].axhline(MC_CUT, ls=":", c="gray")
axs[2].errorbar(xs, tab.b_MLE, yerr=tab.b_err, fmt="o-", c="k", capsize=4); axs[2].set_ylabel(f"b-value (MLE, M≥{MC_CUT})")
axs[2].axhline(1.0, ls=":", c="gray")
for a in axs:
    a.set_xticks(xs); a.set_xticklabels(tab.period, rotation=20); a.grid(alpha=0.2)
plt.tight_layout(); save(fig, "Fig6_temporal_evolution_summary.png"); plt.close(fig)
print("Figures written to", OUT)

# ------------------------------ catalog export ----------------------------- #
cat.to_csv(os.path.join(OUT, f"catalog_M{MC_CUT}_{START_YEAR}plus.csv"), index=False)  # includes along_km
