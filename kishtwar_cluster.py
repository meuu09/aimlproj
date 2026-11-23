"""
Kishtwar-Chamba intermediate-depth cluster: the paper's cluster analysis applied to
our region (USGS ComCat, 1990 to now).

Michailos et al. (2021) analyse "Cluster 1" beneath the Everest region with:
  Fig 3-4  depth cross-sections with depth error bars   -> FigK5 here
  Fig 5    magnitude vs time + cumulative number        -> FigK1
  Fig 6    inter-event times                            -> FigK2
  Fig 7    magnitude-frequency, Mc (max curvature +0.4), b-value  -> FigK3
  Fig 8-9  migration: along/across-strike distance vs time, plus a
           systematic search for sequences within 20 days and 20 km
           of a first event, tested against linear and diffusive
           migration (r = sqrt(4 pi D t), Shapiro et al. 1997)  -> FigK4
  Sec 3.2.1  Bath's-law check (largest minus 2nd-largest magnitude)

Our equivalent of Cluster 1 is the band of deep events between Kishtwar and
Chamba (J&K / Himachal Pradesh). Depth rules (same as him_ides_himachal.py):
  * only events whose depth USGS solved for ("from location" / "constrained by
    depth phases") are used at all; fixed and unflagged depths are dropped
  * intermediate  = depth - depthError >= 40 km (deep even at the top of its error bar)
  * shallow       = depth + depthError <  40 km
  * uncertain     = error bar crosses 40 km, or no depthError given

Run:  .venv/bin/python kishtwar_cluster.py                   (results/excluding_fixed/kishtwar/)
      .venv/bin/python kishtwar_cluster.py --include-fixed   (results/including_fixed/kishtwar/)
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.path import Path

# ----------------------------- PARAMETERS ---------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "data", "usgs_himachal_kashmir_M3plus_1970-2026.csv")
DEPTH_TYPES = os.path.join(HERE, "data", "usgs_depth_type.csv")
# Two cases, chosen on the command line:
#   (default)          excluding_fixed: only USGS-solved depths; fixed/unflagged events dropped
#   --include-fixed    including_fixed: every event, depths as USGS reports them
# Both use the same intermediate rule (depth - depthError >= 40 km).
INCLUDE_FIXED = "--include-fixed" in sys.argv
CASE = "including_fixed" if INCLUDE_FIXED else "excluding_fixed"
USED = "all depths" if INCLUDE_FIXED else "solved depths"   # wording for labels
CASE_LABEL = ("Case: including fixed/unflagged depths" if INCLUDE_FIXED
              else "Case: excluding fixed depths (USGS-solved only)")
OUT = os.path.join(HERE, "results", CASE, "kishtwar")
START_YEAR = 1990
MC_CUT = 4.4                     # regional completeness (authors' method gives 4.0 + 0.4)
DEEP_MIN, DEEP_MAX = 40, 150     # intermediate depth (km)
SOLVED_DEPTH = ("from location", "constrained by depth phases")
# Cluster outline (lon, lat), drawn around the NW-SE Kishtwar-Chamba band of solved
# intermediate-depth events; a few isolated events outside the band are left out.
CLUSTER_POLY = [(75.50, 33.60), (76.10, 33.60), (76.50, 32.95),
                (76.50, 32.75), (75.90, 32.75), (75.50, 33.20)]
SEQ_DAYS, SEQ_KM, SEQ_MIN_EVENTS = 20, 20, 3     # paper's migration-search window
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    """Save a figure into this case's folder, stamped with the case name."""
    fig.text(0.005, 0.995, CASE_LABEL, ha="left", va="top", fontsize=8, color="0.4")
    fig.savefig(os.path.join(OUT, name), dpi=200)
matplotlib.rc("font", size=12)
C_DEEP, C_SHAL, C_REST, C_FIX = "#d62728", "#1f77b4", "#7f7f7f", "#bdbdbd"


# ----------------------------- HELPERS ------------------------------------- #
def max_curv_eqcorrscan(mags, bin_size=0.5):
    """Authors' Mc: eqcorrscan 0.4.2 calc_max_curv (peak |2nd derivative| of 0.5 bins)."""
    min_bin, max_bin = int(min(mags)), int(max(mags) + 1)
    bins = np.arange(min_bin, max_bin + bin_size, bin_size)
    df, bins = np.histogram(mags, bins)
    curvature = np.diff(np.diff(df) / bin_size) / bin_size
    return float(bins[np.argmax(np.abs(curvature))] + bin_size)


def b_lsq(mags, mc):
    """Paper's b-value: straight-line fit to log10(cumulative number) for M >= mc."""
    m = np.sort(np.asarray(mags)[np.asarray(mags) >= mc])
    steps = np.unique(m)
    if len(steps) < 3:
        return np.nan, steps, None
    cum = np.array([(m >= s).sum() for s in steps])
    fit = np.polyfit(steps, np.log10(cum), 1)
    return -fit[0], steps, np.poly1d(fit)


def b_mle(mags, mc, dm=0.1):
    """Aki (1965) / Utsu maximum-likelihood b-value, 1-sigma = b/sqrt(N) (Aki 1965)."""
    m = np.asarray(mags)[np.asarray(mags) >= mc]
    if len(m) < 5:
        return np.nan, np.nan, len(m)
    b = np.log10(np.e) / (m.mean() - (mc - dm / 2))
    return b, b / np.sqrt(len(m)), len(m)


def interevent_days(times):
    t = np.sort(times.values)
    return pd.to_datetime(t[1:]), np.diff(t).astype("timedelta64[s]").astype(float) / 86400.0


def xy_km(lon, lat, lon0, lat0):
    return ((lon - lon0) * 111.32 * np.cos(np.radians(lat0)), (lat - lat0) * 110.57)


# ----------------------------- LOAD ---------------------------------------- #
raw = pd.read_csv(CSV, parse_dates=["time"])
raw["time"] = raw["time"].dt.tz_localize(None)
raw = raw[raw.time.dt.year >= START_YEAR]
raw = raw.merge(pd.read_csv(DEPTH_TYPES)[["id", "depth_type"]], on="id", how="left")
raw["depth_type"] = raw.depth_type.fillna("unknown")
n_before = len(raw)
if not INCLUDE_FIXED:
    raw = raw[raw.depth_type.isin(SOLVED_DEPTH)].copy()   # drop fixed and unflagged depths
raw["in_cluster"] = Path(CLUSTER_POLY).contains_points(raw[["longitude", "latitude"]].values)
raw["intermediate"] = ((raw.depth - raw.depthError) >= DEEP_MIN) & (raw.depth <= DEEP_MAX)
raw["shallow"] = (raw.depth + raw.depthError) < DEEP_MIN
raw["uncertain"] = ~raw.intermediate & ~raw.shallow       # crosses 40 km, or no error
raw = raw.sort_values("time").reset_index(drop=True)

clu_all = raw[raw.in_cluster]                                  # every event used in this case, inside the outline
clu_deep = raw[raw.in_cluster & raw.intermediate]              # our "Cluster 1"
clu_shal = raw[raw.in_cluster & raw.shallow]
clu_unc = raw[raw.in_cluster & raw.uncertain]
rest_deep = raw[~raw.in_cluster & raw.intermediate]            # intermediate elsewhere in box
print(f"{CASE_LABEL} | {START_YEAR}+: {len(raw)} events used, {n_before - len(raw)} fixed/unflagged dropped")
print(f"Cluster outline: used={len(clu_all)}  intermediate={len(clu_deep)}  shallow={len(clu_shal)}  "
      f"uncertain={len(clu_unc)}  | intermediate outside cluster={len(rest_deep)}")
summary = {"N_cluster_used_M3": len(clu_all), "N_cluster_intermediate_M3": len(clu_deep),
           "N_cluster_uncertain_M3": len(clu_unc),
           "N_cluster_intermediate_M>=Mc": int((clu_deep.mag >= MC_CUT).sum()),
           "N_cluster_shallow_M3": len(clu_shal), "N_intermediate_outside": len(rest_deep),
           "share_of_region_intermediate_in_cluster": len(clu_deep) / max(1, len(clu_deep) + len(rest_deep)),
           "depth_mean_km": clu_deep.depth.mean(), "depth_median_km": clu_deep.depth.median(),
           "depth_min_km": clu_deep.depth.min(), "depth_max_km": clu_deep.depth.max(),
           "depthError_median_km": clu_deep.depthError.median()}


def locator(ax):
    poly = np.array(CLUSTER_POLY + [CLUSTER_POLY[0]])
    ax.scatter(raw.longitude, raw.latitude, s=6, c=C_FIX, lw=0)
    ax.scatter(clu_unc.longitude, clu_unc.latitude, s=16, c=C_REST, marker="+")
    ax.scatter(rest_deep.longitude, rest_deep.latitude, s=16, c=C_REST, marker="+")
    ax.scatter(clu_shal.longitude, clu_shal.latitude, s=16, facecolors="none", edgecolors=C_SHAL)
    ax.scatter(clu_deep.longitude, clu_deep.latitude, s=20, c=C_DEEP, marker="X", lw=0)
    ax.plot(poly[:, 0], poly[:, 1], "k-", lw=1)
    ax.set_xlim(74, 79.5); ax.set_ylim(31, 35); ax.set_aspect(1 / np.cos(np.radians(33)))
    ax.yaxis.tick_right(); ax.yaxis.set_label_position("right")
    ax.set_xlabel("Longitude (°)"); ax.set_ylabel("Latitude (°)")


def time_axis(ax):
    ax.set_xlim(pd.Timestamp(f"{START_YEAR}-01-01"), raw.time.max() + pd.Timedelta(days=200))
    ax.xaxis.set_major_locator(mdates.YearLocator(5)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))


# ------------- FIG K1 (paper Fig 5): magnitude vs time + cumulative -------- #
fig = plt.figure(figsize=(13, 6.5))
ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=3, rowspan=2)
for df, col, mk, lab in [(rest_deep, "k", "x", f"Intermediate, outside cluster (N={len(rest_deep)})"),
                         (clu_unc, C_REST, "+", f"Cluster, uncertain: error bar crosses 40 km (N={len(clu_unc)})"),
                         (clu_shal, C_SHAL, "o", f"Cluster, shallow: depth + error < 40 km (N={len(clu_shal)})"),
                         (clu_deep, C_DEEP, "X", f"Cluster, intermediate: depth − error ≥ 40 km (N={len(clu_deep)})")]:
    ax1.vlines(df.time, 0, df.mag, color="gray", lw=0.3, zorder=1)
    ax1.scatter(df.time, df.mag, s=df.mag ** 2 * 3, c="none" if mk == "o" else col,
                edgecolors=None if mk in ("x", "+") else col,
                marker=mk, lw=1.4, alpha=0.85, zorder=3, label=lab)
ax1.axhline(MC_CUT, color="k", ls=":", lw=0.8)
ax1.text(raw.time.min(), MC_CUT + 0.03, f" Mc = {MC_CUT}", fontsize=9)
ax1.set_ylim(3, 5.8); ax1.set_ylabel("Magnitude (mostly mb)"); ax1.set_xlabel("Time (UTC)")
time_axis(ax1); ax1.legend(loc="upper left", fontsize=9)
ax2 = ax1.twinx()
ax2.plot(clu_deep.time, np.arange(1, len(clu_deep) + 1), "k--", lw=2, label="Cluster, intermediate")
ax2.plot(clu_all.time, np.arange(1, len(clu_all) + 1), "k:", lw=2, label=f"Cluster, all events ({USED})")
ax2.plot(rest_deep.time, np.arange(1, len(rest_deep) + 1), "k-", lw=1.2, label="Intermediate, outside")
ax2.set_ylabel("Cumulative number", color="gray"); ax2.tick_params(axis="y", labelcolor="gray")
ax2.legend(loc="lower right", fontsize=9)
locator(plt.subplot2grid((2, 4), (0, 3)))
plt.tight_layout(); save(fig, "FigK1_magnitude_vs_time.png"); plt.close(fig)

# ------------- FIG K2 (paper Fig 6): inter-event times --------------------- #
fig = plt.figure(figsize=(13, 6.5))
ax = plt.subplot2grid((2, 4), (0, 0), colspan=3, rowspan=2)
ie = {}
for lab, df, col, mk in [("cluster_intermediate_M3", clu_deep, C_DEEP, "X"),
                         ("cluster_shallow_M3", clu_shal, C_SHAL, "o"),
                         (f"cluster_intermediate_M{MC_CUT}", clu_deep[clu_deep.mag >= MC_CUT], "k", "s")]:
    t, d = interevent_days(df.time)
    ie[lab] = d
    face = "none" if mk in ("o", "s") else col
    ax.scatter(t, d, s=55, marker=mk, c=face, edgecolors=col, lw=1.4, alpha=0.85,
               label=f"{lab.replace('_', ' ')} (N={len(df)})")
    summary[f"interevent_median_days_{lab}"] = np.median(d) if len(d) else np.nan
    summary[f"interevent_CV_{lab}"] = d.std() / d.mean() if len(d) > 1 else np.nan
    summary[f"interevent_min_days_{lab}"] = d.min() if len(d) else np.nan
ax.set_yscale("log"); ax.set_ylabel("Inter-event time (days)"); ax.set_xlabel("Time (UTC)")
time_axis(ax); ax.legend(loc="lower left", fontsize=9); ax.grid(alpha=0.2, which="both")
locator(plt.subplot2grid((2, 4), (0, 3)))
plt.tight_layout(); save(fig, "FigK2_interevent_times.png"); plt.close(fig)

# ------------- FIG K3 (paper Fig 7): magnitude-frequency + b --------------- #
fig, ax1 = plt.subplots(figsize=(10, 6.5))
bins = np.arange(3.05, raw.mag.max() + 0.2, 0.2)
ax1.hist(raw.mag, bins, histtype="step", fill=True, color="gray", alpha=0.25, lw=1.5, ec="k",
         label=f"Whole region (N={len(raw)})")
ax1.hist(clu_all.mag, bins, histtype="step", lw=1.8, ls="--", ec="k", label=f"Cluster outline, all (N={len(clu_all)})")
ax1.set_xlabel("Magnitude"); ax1.set_ylabel("Frequency (0.2 bins)"); ax1.legend(loc="center left", fontsize=9)
ax2 = ax1.twinx()
bres = []
for lab, df, col, ls, mk in [("Whole region", raw, "k", "-", "+"),
                             (f"Cluster outline, all events ({USED})", clu_all, "#555", "--", "o"),
                             ("Cluster, intermediate", clu_deep, C_DEEP, ":", "X")]:
    m = df.mag.values
    mc_auth = max_curv_eqcorrscan(m) + 0.4
    steps = np.unique(m); cum = np.array([(m >= s).sum() for s in steps])
    ax2.scatter(steps, np.log10(cum), marker=mk, s=30, c="none" if mk == "o" else col, edgecolors=col, lw=1.2)
    b, st, fn = b_lsq(m, MC_CUT)
    bm, bs, n = b_mle(m, MC_CUT)
    if fn is not None:
        ax2.plot(st, fn(st), ls, c=col, lw=2.2,
                 label=f"{lab}: b={b:.2f} (fit), {bm:.2f}±{bs:.2f} (MLE), N={n}")
    bres.append(dict(group=lab, N_total=len(m), Mc_authors_method=mc_auth, Mc_used=MC_CUT,
                     N_above_Mc=n, b_lsq=b, b_mle=bm, b_mle_err=bs))
ax2.axvline(MC_CUT, color="k", lw=0.8, ls=":")
ax2.set_ylabel("log$_{10}$ cumulative number"); ax2.legend(loc="upper right", fontsize=9)
ax1.set_title(f"b-values above M{MC_CUT}; small N, read the MLE ± error", fontsize=11)
plt.tight_layout(); save(fig, "FigK3_magnitude_frequency_bvalue.png"); plt.close(fig)
btab = pd.DataFrame(bres); btab.to_csv(os.path.join(OUT, "kishtwar_bvalues.csv"), index=False)

# ------------- profile axes for the cluster (PCA of epicentres) ------------ #
pts = raw[raw.in_cluster].copy()                       # events used in this case, inside the outline
lon0, lat0 = pts.longitude.mean(), pts.latitude.mean()
x, y = xy_km(pts.longitude, pts.latitude, lon0, lat0)
_, _, vt = np.linalg.svd(np.c_[x - x.mean(), y - y.mean()], full_matrices=False)
strike = vt[0] if vt[0][0] > 0 else -vt[0]             # point the strike axis eastward
across = np.array([-strike[1], strike[0]])             # 90 deg anticlockwise (toward NE)
az = (np.degrees(np.arctan2(strike[0], strike[1])) + 360) % 360
for d in (pts, raw):
    xx, yy = xy_km(d.longitude, d.latitude, lon0, lat0)
    d["along_km"] = np.c_[xx, yy] @ strike
    d["across_km"] = np.c_[xx, yy] @ across
shift = pts.along_km.min()
pts["along_km"] -= shift
summary["profile_azimuth_deg"] = az

# ------------- migration search (paper Sec 3.2.3) -------------------------- #
def hypo_km(a, b):
    xa, ya = xy_km(a.longitude, a.latitude, lon0, lat0)
    xb, yb = xy_km(b.longitude, b.latitude, lon0, lat0)
    return np.sqrt((xa - xb) ** 2 + (ya - yb) ** 2 + (a.depth - b.depth) ** 2)

search = raw[raw.in_cluster & raw.intermediate].reset_index(drop=True)   # paper searched its deep catalog
seqs, members_seen = [], set()
for i, ev in search.iterrows():
    later = search[(search.time > ev.time) & (search.time <= ev.time + pd.Timedelta(days=SEQ_DAYS))]
    later = later[hypo_km(later, ev) <= SEQ_KM]
    ids = tuple([ev.id] + list(later.id))
    if len(ids) >= SEQ_MIN_EVENTS and not set(ids) <= members_seen:
        members_seen |= set(ids)
        grp = search[search.id.isin(ids)].sort_values("time")
        t_s = (grp.time - ev.time).dt.total_seconds().values
        r_m = hypo_km(grp, ev).values * 1000.0
        ok = t_s > 0
        v = np.polyfit(t_s[ok] / 86400, r_m[ok] / 1000, 1)[0] if ok.sum() >= 2 else np.nan
        D = np.max(r_m[ok] ** 2 / (4 * np.pi * t_s[ok])) if ok.any() else np.nan
        seqs.append(dict(start=ev.time, n_events=len(grp), duration_days=t_s.max() / 86400,
                         max_dist_km=r_m.max() / 1000, depths_km=", ".join(f"{z:.0f}" for z in grp.depth),
                         mags=", ".join(f"{m:.1f}" for m in grp.mag),
                         linear_velocity_km_per_day=v, diffusivity_envelope_m2_per_s=D,
                         median_horizontalError_km=grp.horizontalError.median(),
                         median_depthError_km=grp.depthError.median(), ids=";".join(grp.id)))
seq_tab = pd.DataFrame(seqs, columns=["start", "n_events", "duration_days", "max_dist_km", "depths_km", "mags",
                                       "linear_velocity_km_per_day", "diffusivity_envelope_m2_per_s",
                                       "median_horizontalError_km", "median_depthError_km", "ids"])
seq_tab.to_csv(os.path.join(OUT, "kishtwar_sequences.csv"), index=False)
summary["n_sequences_20d_20km"] = len(seq_tab)

# ------------- FIG K4 (paper Figs 8-9): migration plots -------------------- #
fig = plt.figure(figsize=(13, 9))
axm = plt.subplot2grid((3, 3), (0, 0))
axa = plt.subplot2grid((3, 3), (1, 0), colspan=3)
axb = plt.subplot2grid((3, 3), (2, 0), colspan=3, sharex=axa)
deep_mask = pts.intermediate
for ax_, col in [(axa, "across_km"), (axb, "along_km")]:
    ax_.scatter(pts.time[~deep_mask], pts[col][~deep_mask], s=pts.mag[~deep_mask] ** 2 * 4, marker="o",
                c="none", edgecolors="k", lw=0.8, label="shallow or uncertain depth")
    sc = ax_.scatter(pts.time[deep_mask], pts[col][deep_mask], s=pts.mag[deep_mask] ** 2 * 4,
                     c=pts.along_km[deep_mask], cmap="jet", edgecolors="k", lw=0.6, label="intermediate")
    ax_.grid(alpha=0.2)
for s_ in seqs:
    for ax_ in (axa, axb):
        ax_.axvspan(s_["start"], s_["start"] + pd.Timedelta(days=max(s_["duration_days"], 30)), color="orange", alpha=0.25, lw=0)
axa.set_ylabel("Across-strike (km)"); axb.set_ylabel("Along-strike (km)"); axb.set_xlabel("Time (UTC)")
time_axis(axb); axa.legend(loc="upper left", fontsize=9)
plt.setp(axa.get_xticklabels(), visible=False)
axm.scatter(pts.along_km, pts.across_km, s=pts.mag ** 2 * 3, c=np.where(deep_mask, pts.along_km, np.nan),
            cmap="jet", edgecolors="k", lw=0.5)
axm.set_xlabel(f"Along-strike, az {az:.0f}° (km)"); axm.set_ylabel("Across (km)"); axm.set_aspect("equal")
axm.set_title("Map view, profile frame", fontsize=10)
cax = fig.add_axes([0.45, 0.78, 0.4, 0.025])
fig.colorbar(sc, cax=cax, orientation="horizontal", label="Along-strike distance (km), intermediate events")
fig.text(0.45, 0.66, f"Orange bands: sequences with ≥{SEQ_MIN_EVENTS} events within {SEQ_DAYS} days and "
         f"{SEQ_KM} km\n(paper's search window). N sequences = {len(seq_tab)}", fontsize=10)
fig.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.07, hspace=0.35)
save(fig, "FigK4_migration.png"); plt.close(fig)

# ------------- FIG K5 (paper Figs 3-4): depth sections with errors --------- #
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 5.5), gridspec_kw=dict(width_ratios=[1, 1.3]))
for d, col, lab in [(pts[pts.shallow], C_SHAL, "shallow: depth + error < 40 km"),
                    (pts[pts.uncertain], C_REST, "uncertain: crosses 40 km / no error"),
                    (pts[deep_mask], C_DEEP, "intermediate: depth − error ≥ 40 km")]:
    a1.errorbar(d.along_km, d.depth, yerr=d.depthError.fillna(0), fmt="o", ms=5, c=col, ecolor=col,
                alpha=0.8, elinewidth=0.8, capsize=2, label=lab)
    a2.errorbar(d.time, d.depth, yerr=d.depthError.fillna(0), fmt="o", ms=5, c=col, ecolor=col,
                alpha=0.8, elinewidth=0.8, capsize=2, label=lab)
for a in (a1, a2):
    a.invert_yaxis(); a.axhline(DEEP_MIN, ls=":", c="k", lw=0.8); a.set_ylabel("Depth (km)"); a.grid(alpha=0.2)
a1.set_xlabel(f"Along-strike distance, az {az:.0f}° (km)"); a1.legend(fontsize=9, loc="lower left")
a2.set_xlabel("Time (UTC)"); time_axis(a2)
a1.set_title("Depth section along the cluster (bars = USGS depthError)", fontsize=10)
a2.set_title("Depth through time", fontsize=10)
plt.tight_layout(); save(fig, "FigK5_depth_sections.png"); plt.close(fig)

# ------------- depth vs catalog era (is a depth "trend" real?) -------------- #
# Location quality changes over time (networks, methods), so compare depth AND its error.
for lab, d in [("before2010", clu_deep[clu_deep.time.dt.year < 2010]),
               ("2010on", clu_deep[clu_deep.time.dt.year >= 2010])]:
    summary[f"N_intermediate_{lab}"] = len(d)
    summary[f"depth_median_km_{lab}"] = d.depth.median()
    summary[f"depthError_median_km_{lab}"] = d.depthError.median()

# ------------- Bath's law check (paper Sec 3.2.1) ------------------------- #
top = np.sort(clu_deep.mag.values)[::-1]
summary["largest_mag"] = top[0]
summary["largest_minus_second"] = top[0] - top[1] if len(top) > 1 else np.nan
for _, s_ in seq_tab.iterrows():
    mm = sorted(map(float, s_.mags.split(", ")), reverse=True)
    print(f"Sequence {s_.start:%Y-%m-%d}: N={s_.n_events}, mags {s_.mags}, depths {s_.depths_km} km, "
          f"dM(1st-2nd)={mm[0] - mm[1]:.1f}, v={s_.linear_velocity_km_per_day:.2f} km/d, "
          f"D<={s_.diffusivity_envelope_m2_per_s:.1f} m2/s, median errors h/z = "
          f"{s_.median_horizontalError_km:.1f}/{s_.median_depthError_km:.1f} km")

# ------------- outputs ------------------------------------------------------ #
raw[raw.in_cluster].drop(columns=["along_km", "across_km"]).merge(
    pts[["id", "along_km", "across_km"]], on="id", how="left").to_csv(
    os.path.join(OUT, "kishtwar_catalog.csv"), index=False)
pd.Series(summary).to_csv(os.path.join(OUT, "kishtwar_summary.csv"), header=["value"])
print(pd.Series(summary).round(3).to_string())
print(btab.round(2).to_string(index=False))
print("Outputs in", OUT)
