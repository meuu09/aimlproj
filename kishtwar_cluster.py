"""Kishtwar-Chamba intermediate-depth cluster: the paper's cluster analysis applied to our region."""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.path import Path

from seis_utils import (CLUSTER_POLY, DEEP_MIN, b_lsq, b_mle, case_label, case_name, cv, interevent_days,
                        load_catalog, max_curv_eqcorrscan, parse_args, results_base, xy_km)

ARGS = parse_args("Kishtwar-Chamba intermediate-depth cluster analysis")
CATALOG, INCLUDE_FIXED = ARGS.catalog, ARGS.include_fixed
CASE = case_name(ARGS)
USED = "all depths" if INCLUDE_FIXED else "solved depths"
CASE_LABEL = case_label(ARGS)
OUT = os.path.join(results_base(ARGS), CASE, "kishtwar")
START_YEAR = 1990
MC_CUT = 4.4
MAG_WORD = "mostly mb" if CATALOG == "usgs" else "mixed types: mb, ML, MD, Mw"
SEQ_DAYS, SEQ_KM, SEQ_MIN_EVENTS = 20, 20, 3
os.makedirs(OUT, exist_ok=True)


def save(fig, name):
    """Save a figure into this case's folder, stamped with the case name."""
    fig.text(0.005, 0.995, CASE_LABEL, ha="left", va="top", fontsize=8, color="0.4")
    fig.savefig(os.path.join(OUT, name), dpi=200)


matplotlib.rc("font", size=12)
C_DEEP, C_SHAL, C_REST, C_FIX = "#d62728", "#1f77b4", "#7f7f7f", "#bdbdbd"

raw, load_info = load_catalog(CATALOG, INCLUDE_FIXED, START_YEAR, ARGS.drop_duplicates)
n_before = load_info["after_duplicates"]
raw["in_cluster"] = Path(CLUSTER_POLY).contains_points(raw[["longitude", "latitude"]].values)
raw["intermediate"] = raw.depth_class == "intermediate"
raw["shallow"] = raw.depth_class == "shallow"
raw["deep"] = raw.depth_class == "deep"
raw["uncertain"] = raw.depth_class == "uncertain"

clu_all = raw[raw.in_cluster]
clu_deep = raw[raw.in_cluster & raw.intermediate]
clu_shal = raw[raw.in_cluster & raw.shallow]
clu_unc = raw[raw.in_cluster & raw.uncertain]
rest_deep = raw[~raw.in_cluster & raw.intermediate]
print(f"{CASE_LABEL} | {START_YEAR}+: {len(raw)} events used, {n_before - len(raw)} fixed/unflagged dropped")
print(f"Cluster outline: used={len(clu_all)}  intermediate={len(clu_deep)}  shallow={len(clu_shal)}  "
      f"uncertain={len(clu_unc)}  | intermediate outside cluster={len(rest_deep)}")
summary = {"N_cluster_used_M3": len(clu_all), "N_cluster_intermediate_M3": len(clu_deep),
           "N_cluster_uncertain_M3": len(clu_unc),
           "N_cluster_intermediate_M>=Mc": int((clu_deep.mag >= MC_CUT).sum()),
           "N_cluster_shallow_M3": len(clu_shal),
           "N_cluster_deeper_than_150km": int(raw[raw.in_cluster & raw.deep].shape[0]),
           "lowest_magnitude_used": raw.mag.min(), "N_intermediate_outside": len(rest_deep),
           "share_of_region_intermediate_in_cluster": len(clu_deep) / max(1, len(clu_deep) + len(rest_deep)),
           "depth_mean_km": clu_deep.depth.mean(), "depth_median_km": clu_deep.depth.median(),
           "depth_min_km": clu_deep.depth.min(), "depth_max_km": clu_deep.depth.max(),
           "depthError_median_km": clu_deep.depthError.median()}


def locator(ax):
    """Small map of the region with the cluster outline."""
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
    """Time axis from START_YEAR to just after the last event, ticks every 5 years."""
    ax.set_xlim(pd.Timestamp(f"{START_YEAR}-01-01"), raw.time.max() + pd.Timedelta(days=200))
    ax.xaxis.set_major_locator(mdates.YearLocator(5)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))


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
ax1.set_ylim(min(3.0, np.floor(raw.mag.min() * 10) / 10 - 0.1), max(5.8, raw.mag.max() + 0.2))
ax1.set_ylabel(f"Magnitude ({MAG_WORD})"); ax1.set_xlabel("Time (UTC)")
time_axis(ax1); ax1.legend(loc="upper left", fontsize=9)
ax2 = ax1.twinx()
ax2.plot(clu_deep.time, np.arange(1, len(clu_deep) + 1), "k--", lw=2, label="Cluster, intermediate")
ax2.plot(clu_all.time, np.arange(1, len(clu_all) + 1), "k:", lw=2, label=f"Cluster, all events ({USED})")
ax2.plot(rest_deep.time, np.arange(1, len(rest_deep) + 1), "k-", lw=1.2, label="Intermediate, outside")
ax2.set_ylabel("Cumulative number", color="gray"); ax2.tick_params(axis="y", labelcolor="gray")
ax2.legend(loc="lower right", fontsize=9)
locator(plt.subplot2grid((2, 4), (0, 3)))
plt.tight_layout(); save(fig, "FigK1_magnitude_vs_time.png"); plt.close(fig)

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
    summary[f"interevent_CV_{lab}"] = cv(d)
    summary[f"interevent_min_days_{lab}"] = d.min() if len(d) else np.nan
ax.set_yscale("log"); ax.set_ylabel("Inter-event time (days)"); ax.set_xlabel("Time (UTC)")
time_axis(ax); ax.legend(loc="lower left", fontsize=9); ax.grid(alpha=0.2, which="both")
locator(plt.subplot2grid((2, 4), (0, 3)))
plt.tight_layout(); save(fig, "FigK2_interevent_times.png"); plt.close(fig)

fig, ax1 = plt.subplots(figsize=(10, 6.5))
bins = np.arange(np.floor(raw.mag.min() * 10) / 10 - 0.05, raw.mag.max() + 0.2, 0.2)
ax1.hist(raw.mag, bins, histtype="step", fill=True, color="gray", alpha=0.25, lw=1.5, ec="k",
         label=f"Whole region (N={len(raw)})")
ax1.hist(clu_all.mag, bins, histtype="step", lw=1.8, ls="--", ec="k", label=f"Cluster outline, all (N={len(clu_all)})")
ax1.set_xlabel(f"Magnitude ({MAG_WORD})"); ax1.set_ylabel("Frequency (0.2 bins)"); ax1.legend(loc="center left", fontsize=9)
ax2 = ax1.twinx()
bres = []
for lab, df, col, ls, mk in [("Whole region", raw, "k", "-", "+"),
                             (f"Cluster outline, all events ({USED})", clu_all, "#555", "--", "o"),
                             ("Cluster, intermediate", clu_deep, C_DEEP, ":", "X")]:
    m = df.mag.values
    if len(m) == 0:
        bres.append(dict(group=lab, N_total=0)); continue
    mc_auth = max_curv_eqcorrscan(m) + 0.4
    steps = np.unique(m); cum = np.array([(m >= s).sum() for s in steps])
    ax2.scatter(steps, np.log10(cum), marker=mk, s=30, c="none" if mk == "o" else col, edgecolors=col, lw=1.2)
    b, st, _, fn = b_lsq(m, MC_CUT)
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

pts = raw[raw.in_cluster].copy()
lon0, lat0 = pts.longitude.mean(), pts.latitude.mean()
x, y = xy_km(pts.longitude, pts.latitude, lon0, lat0)
_, _, vt = np.linalg.svd(np.c_[x - x.mean(), y - y.mean()], full_matrices=False)
strike = vt[0] if vt[0][0] > 0 else -vt[0]
across = np.array([-strike[1], strike[0]])
az = (np.degrees(np.arctan2(strike[0], strike[1])) + 360) % 360
pts["along_km"] = np.c_[x, y] @ strike
pts["across_km"] = np.c_[x, y] @ across
pts["along_km"] -= pts.along_km.min()
summary["profile_azimuth_deg"] = az


def hypo_km(a, b):
    """3-D hypocentral distance (km) between events a and b."""
    xa, ya = xy_km(a.longitude, a.latitude, lon0, lat0)
    xb, yb = xy_km(b.longitude, b.latitude, lon0, lat0)
    return np.sqrt((xa - xb) ** 2 + (ya - yb) ** 2 + (a.depth - b.depth) ** 2)


search = clu_deep.reset_index(drop=True)
seqs, members_seen = [], set()
for _, ev in search.iterrows():
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
                         linear_velocity_km_per_day=v, diffusivity_min_envelope_m2_per_s=D,
                         median_horizontalError_km=grp.horizontalError.median(),
                         median_depthError_km=grp.depthError.median(), ids=";".join(grp.id)))
seq_tab = pd.DataFrame(seqs, columns=["start", "n_events", "duration_days", "max_dist_km", "depths_km", "mags",
                                       "linear_velocity_km_per_day", "diffusivity_min_envelope_m2_per_s",
                                       "median_horizontalError_km", "median_depthError_km", "ids"])
seq_tab.to_csv(os.path.join(OUT, "kishtwar_sequences.csv"), index=False)
summary["n_sequences_20d_20km"] = len(seq_tab)

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

for lab, d in [("before2010", clu_deep[clu_deep.time.dt.year < 2010]),
               ("2010on", clu_deep[clu_deep.time.dt.year >= 2010])]:
    summary[f"N_intermediate_{lab}"] = len(d)
    summary[f"depth_median_km_{lab}"] = d.depth.median()
    summary[f"depthError_median_km_{lab}"] = d.depthError.median()

top = np.sort(clu_deep.mag.values)[::-1]
summary["largest_mag"] = top[0] if len(top) else np.nan
summary["largest_minus_second"] = top[0] - top[1] if len(top) > 1 else np.nan
for _, s_ in seq_tab.iterrows():
    mm = sorted(map(float, s_.mags.split(", ")), reverse=True)
    print(f"Sequence {s_.start:%Y-%m-%d}: N={s_.n_events}, mags {s_.mags}, depths {s_.depths_km} km, "
          f"dM(1st-2nd)={mm[0] - mm[1]:.1f}, v={s_.linear_velocity_km_per_day:.2f} km/d, "
          f"D>={s_.diffusivity_min_envelope_m2_per_s:.1f} m2/s, median errors h/z = "
          f"{s_.median_horizontalError_km:.1f}/{s_.median_depthError_km:.1f} km")

pts.to_csv(os.path.join(OUT, "kishtwar_catalog.csv"), index=False)
pd.Series(summary).to_csv(os.path.join(OUT, "kishtwar_summary.csv"), header=["value"])
print(pd.Series(summary).round(3).to_string())
print(btab.round(2).to_string(index=False))
print("Outputs in", OUT)
