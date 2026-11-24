"""
Spatial earthquake clusters, found separately for each depth group, over the
complete dataset (all years in the catalog).

Depth groups (same rule as the other scripts):
  intermediate  depth - depthError >= 40 km
  shallow       depth + depthError <  40 km
  uncertain     error bar crosses 40 km, or no depth error
Clusters are found with DBSCAN (Ester et al., 1996), an unsupervised machine-learning
method: an event belongs to a cluster if at least MIN_EVENTS events (itself included)
lie within EPS_KM of it on the map, and clusters grow through such dense neighbours.
Events in no dense group are "noise". The number of clusters is not set in advance.
EPS_KM = 30 km is about the size of USGS horizontal location errors in this region.

For each cluster we report size, centre, extent, depth range and error, time span,
magnitude range, and (when enough events) the b-value and inter-event CV.

Run (each combination writes its own folder):
  .venv/bin/python depth_clusters.py                          USGS, excluding fixed
  .venv/bin/python depth_clusters.py --include-fixed          USGS, including fixed
  .venv/bin/python depth_clusters.py --catalog combined [...]  USGS + ISC Bulletin
  add --drop-duplicates (combined only) to leave out rows flagged as possible duplicates
Windows (PowerShell): .venv\\Scripts\\python depth_clusters.py [same options]

A fourth group, "deep" (depth > 150 km), is kept apart; before it was counted as uncertain.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

from seis_utils import (BOX, DEEP_MIN, SOLVED_DEPTH, b_mle, case_label, case_name, load_catalog, parse_args,
                        results_base)
from seis_utils import xy_km as km_xy

# ----------------------------- PARAMETERS ---------------------------------- #
ARGS = parse_args("DBSCAN clusters per depth group")
INCLUDE_FIXED, CATALOG = ARGS.include_fixed, ARGS.catalog
OUT = os.path.join(results_base(ARGS), case_name(ARGS), "depth_clusters")
CASE_LABEL = case_label(ARGS) + " | all years"
EPS_KM, MIN_EVENTS = 30.0, 4          # DBSCAN neighbourhood radius and minimum cluster size
MC_CUT = 4.4                          # completeness used for b-values (USGS mb, authors' method)
PROFILE_AZ = 127.0                    # arc-parallel azimuth (PCA of the USGS M>=4.4 epicentres, him_ides_himachal.py)
GROUPS = [("intermediate", "#c62828"), ("uncertain", "#757575"), ("shallow", "#1565c0"), ("deep", "#000000")]
os.makedirs(OUT, exist_ok=True)
matplotlib.rc("font", size=11)


def load():
    c, info = load_catalog(CATALOG, INCLUDE_FIXED, None, ARGS.drop_duplicates)
    return c.rename(columns={"depth_class": "group"}), info["in_box_with_mag"]


cat, n_all = load()
print(f"{CASE_LABEL}: {n_all} events in box, {len(cat)} used, "
      f"{(cat.time.min()):%Y}-{(cat.time.max()):%Y}")
print(cat.group.value_counts().to_string())

rows = []
cat["cluster"] = -1
for g, _ in GROUPS:
    sub = cat[cat.group == g]
    if len(sub) < MIN_EVENTS:
        continue
    X = np.radians(sub[["latitude", "longitude"]].values)
    lab = DBSCAN(eps=EPS_KM / 6371.0, min_samples=MIN_EVENTS, metric="haversine").fit_predict(X)
    # number clusters by size, largest = 1
    order = pd.Series(lab[lab >= 0]).value_counts().index.tolist()
    remap = {old: i + 1 for i, old in enumerate(order)}
    cat.loc[sub.index, "cluster"] = [remap.get(l, -1) for l in lab]
    for old in order:
        k = remap[old]
        cl = sub[lab == old]
        lat0, lon0 = cl.latitude.mean(), cl.longitude.mean()
        x, y = km_xy(cl.longitude, cl.latitude, lon0, lat0)
        span = float(np.max(np.hypot(x.values[:, None] - x.values[None, :], y.values[:, None] - y.values[None, :])))
        gaps = np.diff(cl.time.values).astype("timedelta64[s]").astype(float) / 86400
        bm, bs, nb = b_mle(cl.mag.values, MC_CUT)
        rows.append(dict(depth_group=g, cluster=k, n_events=len(cl), n_M_ge_Mc=int((cl.mag >= MC_CUT).sum()),
                         lat_centre=round(lat0, 3), lon_centre=round(lon0, 3), max_extent_km=round(span, 1),
                         depth_median_km=round(cl.depth.median(), 1), depth_min_km=cl.depth.min(),
                         depth_max_km=cl.depth.max(), depthError_median_km=round(cl.depthError.median(), 1),
                         first=cl.time.min().date(), last=cl.time.max().date(),
                         mag_min=cl.mag.min(), mag_max=cl.mag.max(),
                         b_mle=round(bm, 2) if bm == bm else None, b_err=round(bs, 2) if bs == bs else None,
                         interevent_CV=round(gaps.std() / gaps.mean(), 2) if len(gaps) > 2 else None,
                         share_fixed_or_unflagged=round((~cl.depth_type.isin(SOLVED_DEPTH)).mean(), 2),
                         n_possible_duplicates=int(cl.possible_duplicate_of.notna().sum()),
                         sources=";".join(sorted(cl.source.astype(str).unique()))))
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, "depth_clusters.csv"), index=False)
cat.to_csv(os.path.join(OUT, "events_with_clusters.csv"), index=False)

# ----------------------------- MAP FIGURE ---------------------------------- #
fig, axs = plt.subplots(1, 4, figsize=(22, 6.2), sharex=True, sharey=True)
pal = plt.get_cmap("tab10")
for ax, (g, col) in zip(axs, GROUPS):
    sub = cat[cat.group == g]
    ax.scatter(cat.longitude, cat.latitude, s=3, c="0.85", lw=0)
    noise = sub[sub.cluster < 0]
    ax.scatter(noise.longitude, noise.latitude, s=14, facecolors="none", edgecolors="0.45", lw=0.7,
               label=f"not in a cluster ({len(noise)})")
    for k in sorted(sub.cluster[sub.cluster > 0].unique()):
        cl = sub[sub.cluster == k]
        ax.scatter(cl.longitude, cl.latitude, s=10 + cl.mag ** 2 * 4, color=pal((k - 1) % 10),
                   edgecolors="k", lw=0.4, alpha=0.9, label=f"cluster {k} ({len(cl)})")
        ax.annotate(str(k), (cl.longitude.mean(), cl.latitude.mean()), fontsize=11, weight="bold",
                    ha="center", va="center", color="k",
                    bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="k", lw=0.6))
    nclu = sub.cluster[sub.cluster > 0].nunique()
    title = {"intermediate": "Intermediate: depth − error ≥ 40 km",
             "uncertain": "Uncertain: error bar crosses 40 km",
             "shallow": "Shallow: depth + error < 40 km",
             "deep": "Deeper than 150 km"}[g]
    ax.set_title(f"{title}\n{len(sub)} events, {nclu} clusters", fontsize=11, color=col)
    ax.set_xlim(BOX["lon"][0] - 0.1, BOX["lon"][1] + 0.1); ax.set_ylim(BOX["lat"][0] - 0.1, BOX["lat"][1] + 0.1)
    ax.set_aspect(1 / np.cos(np.radians(33))); ax.set_xlabel("Longitude (°)"); ax.grid(alpha=0.2)
    if len(sub):
        ax.legend(fontsize=7.5, loc="upper right", ncol=1 if nclu < 8 else 2, framealpha=0.85)
axs[0].set_ylabel("Latitude (°)")
fig.suptitle(f"DBSCAN clusters by depth group (radius {EPS_KM:.0f} km, ≥{MIN_EVENTS} events)   ·   {CASE_LABEL}",
             fontsize=11)
plt.tight_layout(); fig.savefig(os.path.join(OUT, "DepthClusters_map.png"), dpi=200); plt.close(fig)

# --------------------- DEPTH VIEW: clusters in depth ----------------------- #
fig, ax = plt.subplots(figsize=(13, 5.5))
lat0, lon0 = cat.latitude.mean(), cat.longitude.mean()
x, y = km_xy(cat.longitude, cat.latitude, lon0, lat0)
az = np.radians(PROFILE_AZ)                             # regional arc-parallel profile
cat["along_km"] = x * np.sin(az) + y * np.cos(az)
for g, col in GROUPS:
    sub = cat[cat.group == g]
    inc = sub.cluster > 0
    ax.errorbar(sub.along_km[inc], sub.depth[inc], yerr=sub.depthError[inc].fillna(0), fmt="o", ms=5,
                color=col, ecolor=col, elinewidth=0.7, alpha=0.85, label=f"{g}, in a cluster")
    ax.scatter(sub.along_km[~inc], sub.depth[~inc], s=12, facecolors="none", edgecolors=col, lw=0.6, alpha=0.6,
               label=f"{g}, not clustered")
ax.axhline(DEEP_MIN, ls=":", c="k", lw=0.8)
ax.invert_yaxis(); ax.set_ylabel("Depth (km)"); ax.set_xlabel(f"Distance along the arc, azimuth {PROFILE_AZ:.0f}° (km)")
ax.legend(fontsize=8, ncol=3, loc="lower left"); ax.grid(alpha=0.2)
ax.set_title(f"All events projected onto one arc-parallel line; bars = depth error   ·   {CASE_LABEL}", fontsize=10)
plt.tight_layout(); fig.savefig(os.path.join(OUT, "DepthClusters_section.png"), dpi=200); plt.close(fig)

pd.set_option("display.width", 250)
print(tab.drop(columns=["sources"]).to_string(index=False) if len(tab) else "no clusters")
print("Outputs in", OUT)
