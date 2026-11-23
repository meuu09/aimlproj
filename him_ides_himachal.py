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
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# ----------------------------- PARAMETERS ---------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "data", "usgs_himachal_kashmir_M3plus_1970-2026.csv")
# Two cases, chosen on the command line:
#   (default)          excluding_fixed: only USGS-solved depths; fixed/unflagged events dropped
#   --include-fixed    including_fixed: every event, depths as USGS reports them
# Both use the same intermediate rule (depth - depthError >= 40 km).
# Catalog, chosen on the command line: "--catalog usgs" (default) or "--catalog combined"
# (USGS + ISC Bulletin, data/combined_usgs_isc_catalog.csv from build_combined_catalog.py).
CATALOG = sys.argv[sys.argv.index("--catalog") + 1] if "--catalog" in sys.argv else "usgs"
COMBINED = os.path.join(HERE, "data", "combined_usgs_isc_catalog.csv")
INCLUDE_FIXED = "--include-fixed" in sys.argv
CASE = "including_fixed" if INCLUDE_FIXED else "excluding_fixed"
USED = "all depths" if INCLUDE_FIXED else "solved depths"   # wording for labels
CASE_LABEL = (("USGS + ISC | " if CATALOG == "combined" else "USGS | ")
              + ("Case: including fixed/unflagged depths" if INCLUDE_FIXED
                 else "Case: excluding fixed depths (solved only)"))
OUT = os.path.join(HERE, "results", CASE) if CATALOG == "usgs" else os.path.join(HERE, "results", "combined_usgs_isc", CASE)
START_YEAR = 1990          # analysis period
MC_CUT = 4.4               # magnitude cutoff used for rate / inter-event / b-value
DEEP_MIN, DEEP_MAX = 40, 150   # "intermediate depth" range (km)
DEPTH_TYPES = os.path.join(HERE, "data", "usgs_depth_type.csv")
# USGS depth-type values that mean the depth was computed. Only these events are used;
# "operator assigned" (fixed by an analyst) and unflagged events are dropped.
SOLVED_DEPTH = ("from location", "constrained by depth phases")
BOX = dict(lat=(31.0, 35.0), lon=(74.0, 79.5))   # region the catalog was cut to
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    """Save a figure into this case's folder, stamped with the case name."""
    fig.text(0.005, 0.995, CASE_LABEL, ha="left", va="top", fontsize=8, color="0.4")
    fig.savefig(os.path.join(OUT, name), dpi=200)

C_DEEP, C_SHAL, C_FIX = "#d62728", "#1f77b4", "#9a9a9a"
matplotlib.rc("font", size=13)

# ----------------------------- LOAD / PREPARE ------------------------------ #
if CATALOG == "combined":
    raw = pd.read_csv(COMBINED)                                  # depth_type already included
    raw["time"] = pd.to_datetime(raw.time, utc=True, format="ISO8601").dt.tz_localize(None)
    raw = raw[raw.time.dt.year >= START_YEAR].copy()
else:
    raw = pd.read_csv(CSV, parse_dates=["time"])
    raw["time"] = raw["time"].dt.tz_localize(None)
    raw = raw[raw.time.dt.year >= START_YEAR].copy()
    raw = raw.merge(pd.read_csv(DEPTH_TYPES)[["id", "depth_type"]], on="id", how="left")
if "horizontalError" not in raw:
    raw["horizontalError"] = np.nan
raw["depth_type"] = raw.depth_type.fillna("unknown")
n_before = len(raw)
if not INCLUDE_FIXED:
    raw = raw[raw.depth_type.isin(SOLVED_DEPTH)].copy()   # drop fixed and unflagged depths
print(f"{CASE_LABEL} | {START_YEAR}+: {n_before} events, {len(raw)} used, "
      f"{n_before - len(raw)} fixed/unflagged dropped")
top, bottom = raw.depth - raw.depthError, raw.depth + raw.depthError   # NaN error -> uncertain
raw["depth_class"] = np.select(
    [(top >= DEEP_MIN) & (raw.depth <= DEEP_MAX), bottom < DEEP_MIN],
    ["intermediate (depth - error >= 40 km)", "shallow (depth + error < 40 km)"],
    default="uncertain (error bar crosses 40 km, or no error)")
cat = raw[raw.mag >= MC_CUT].sort_values("time").reset_index(drop=True)
deep = cat[cat.depth_class.str.startswith("intermediate")]
shal = cat[cat.depth_class.str.startswith("shallow")]
unc = cat[cat.depth_class.str.startswith("uncertain")]
print(f"Catalog {START_YEAR}+, M>={MC_CUT}: N={len(cat)}  intermediate={len(deep)}  "
      f"shallow={len(shal)}  uncertain={len(unc)}")
print(cat.depth_class.value_counts().to_string())


def locator(ax, df_groups):
    for df, col, mk in df_groups:
        ax.scatter(df.longitude, df.latitude, s=14, c=col, marker=mk, alpha=0.8, lw=0.6)
    ax.set_xlim(BOX["lon"][0] - 0.3, BOX["lon"][1] + 0.3)
    ax.set_ylim(BOX["lat"][0] - 0.3, BOX["lat"][1] + 0.3)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.yaxis.set_label_position("right"); ax.yaxis.tick_right()
    ax.tick_params(top=True, left=True)


def max_curv(mags, step=0.1):
    """Mc by maximum curvature, Wiemer & Wyss (2000): the most frequent magnitude bin."""
    mags = np.asarray(mags)
    bins = np.arange(np.floor(mags.min()) - 1, np.ceil(mags.max()) + 1, step)
    counts, edges = np.histogram(mags, bins=bins)
    return round(float(edges[np.argmax(counts)]), 2)


def max_curv_eqcorrscan(mags, bin_size=0.5):
    """The authors' Mc: eqcorrscan 0.4.2 calc_max_curv, copied verbatim
    (peak |second derivative| of the 0.5-bin histogram). The paper adds +0.4 to it."""
    min_bin, max_bin = int(min(mags)), int(max(mags) + 1)
    bins = np.arange(min_bin, max_bin + bin_size, bin_size)
    df, bins = np.histogram(mags, bins)
    curvature = np.diff(np.diff(df) / bin_size) / bin_size
    return float(bins[np.argmax(np.abs(curvature))] + bin_size)


def b_lsq(mags, mc):
    """Original method: straight-line fit to log10(cumulative number) above mc."""
    m = np.sort(np.asarray(mags)[np.asarray(mags) >= mc])
    steps = np.unique(m)
    cum = np.array([(m >= s).sum() for s in steps])
    fit = np.polyfit(steps, np.log10(cum), 1)
    return -fit[0], steps, cum, np.poly1d(fit)


def b_mle(mags, mc, dm=0.1):
    """Aki-Utsu maximum-likelihood b-value (with binning correction) and 1-sigma."""
    m = np.asarray(mags)[np.asarray(mags) >= mc]
    if len(m) < 5:
        return np.nan, np.nan, len(m)
    b = np.log10(np.e) / (m.mean() - (mc - dm / 2))
    return b, b / np.sqrt(len(m)), len(m)


# ----------------- FIG A: magnitude vs time + cumulative ------------------- #
fig = plt.figure(figsize=(13, 6.5))
ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=3, rowspan=2)
for df, col, mk, lab in [(unc, C_FIX, "+", f"Uncertain: error bar crosses 40 km (N={len(unc)})"),
                         (shal, C_SHAL, "o", f"Shallow: depth + error < 40 km (N={len(shal)})"),
                         (deep, C_DEEP, "X", f"Intermediate: depth − error ≥ 40 km (N={len(deep)})")]:
    ax1.vlines(df.time, 0, df.mag, color="gray", lw=0.3, zorder=1)
    ax1.scatter(df.time, df.mag, s=df.mag**2 * 2.2, c=col if mk != "o" else "none",
                edgecolors=col, marker=mk, alpha=0.8, lw=1.4, zorder=3, label=lab)
ax1.set_ylim(MC_CUT - 0.2, cat.mag.max() + 0.3)
ax1.set_xlim(pd.Timestamp(f"{START_YEAR}-01-01"), cat.time.max() + pd.Timedelta(days=200))
ax1.set_ylabel("Magnitude"); ax1.set_xlabel("Time (UTC)")
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
    t = df.time.values
    return pd.to_datetime(t[1:]), np.diff(t).astype("timedelta64[s]").astype(float) / 86400.0

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
      f"min={da.min():.3f} d, max={da.max():.0f} d | CV={da.std()/da.mean():.2f}")
print(f"Inter-event time, intermediate: median={np.median(dd):.1f} d, mean={dd.mean():.1f} d")

# ---------------------- FIG C: magnitude-frequency + b --------------------- #
mags_all = raw.mag.values                      # solved-depth M>=3 catalog (1990+) to locate Mc
mc_curv = max_curv(mags_all)
fig, ax1 = plt.subplots(figsize=(9.5, 6.5))
bins = np.arange(3.05, raw.mag.max() + 0.2, 0.2)
ax1.hist(raw.mag, bins, histtype="step", fill=True, color="gray", alpha=0.25, lw=2, ec="k",
         label=f"All events ({USED}), M≥3 since {START_YEAR} (N={len(raw)})")
ax1.set_ylabel("Frequency"); ax1.set_xlabel("Magnitude")
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
    if len(steps) >= 2:
        ax2.plot(steps, fn(steps), ls, c=col, lw=2.4,
                 label=f"{lab}: b={b:.2f} (fit), {bm:.2f}±{bs:.2f} (MLE)\nM≥{MC_CUT}, N={n}")
ax2.axvline(MC_CUT, color="k", lw=0.8, ls=":")
ax2.set_ylabel("log$_{10}$ cumulative number")
ax2.legend(loc="upper right", fontsize=10)
ax1.set_title(f"Max-curvature Mc of this catalog = {mc_curv}; fit threshold M≥{MC_CUT}", fontsize=11)
plt.tight_layout(); save(fig, "Fig3_magnitude_frequency_bvalue.png"); plt.close(fig)
print(f"Max-curvature Mc (M>=3, {START_YEAR}+): Wiemer-Wyss mode = {mc_curv}; "
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
cb = plt.colorbar(sc, ax=ax); cb.set_label("Solved depth (km)")
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
ax.set_title(f"USGS-solved depths only, M≥{MC_CUT} (N={len(meas)}); fixed/unflagged depths not used", fontsize=11)
fig.colorbar(p, ax=ax, shrink=0.6, label="Depth (km)"); ax.view_init(elev=22, azim=-65)
plt.tight_layout(); save(fig, "Fig5_3D_hypocentres.png"); plt.close(fig)

# ------------- FIG F: temporal evolution summary (rate, Mc, b) ------------- #
periods = [(1990, 1999), (2000, 2009), (2010, 2019), (2020, 2026)]
rows = []
for a, b_ in periods:
    seg = raw[(raw.time.dt.year >= a) & (raw.time.dt.year <= b_)]
    yrs = min(b_, 2026) - a + 1 if b_ < 2026 else (2026 - a) + 10.5 / 12
    n44 = (seg.mag >= MC_CUT).sum()
    bm, bs, n = b_mle(seg.mag.values, MC_CUT)
    rows.append(dict(period=f"{a}-{b_}", Mc_maxcurv=max_curv(seg.mag.values), N_M3=len(seg), N_Mc=n44,
                     rate_per_yr=n44 / yrs, b_MLE=bm, b_err=bs))
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
