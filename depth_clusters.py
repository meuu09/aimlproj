"""Spatial DBSCAN clusters by depth-quality class, using an explicit M floor.

DBSCAN uses only horizontal coordinates and can chain nearby neighborhoods into
large regions. Results are descriptive density groups, not time-dependent or
uncertainty-aware physical clusters. Cluster b-values are omitted when a cluster
mixes magnitude types or sources.

Run from Windows PowerShell in the repository root:
    .\\.venv\\Scripts\\python.exe depth_clusters.py
    .\\.venv\\Scripts\\python.exe depth_clusters.py --catalog combined-v2
    .\\.venv\\Scripts\\python.exe depth_clusters.py --catalog combined-v2 --min-mag 0
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

# ----------------------------- PARAMETERS ---------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE_FIXED = "--include-fixed" in sys.argv
CATALOG = sys.argv[sys.argv.index("--catalog") + 1] if "--catalog" in sys.argv else "usgs"
MIN_MAG = float(sys.argv[sys.argv.index("--min-mag") + 1]) if "--min-mag" in sys.argv else 3.0
CASE = "including_fixed" if INCLUDE_FIXED else "excluding_fixed"
if CATALOG == "usgs":
    BASE = os.path.join(HERE, "results")
elif CATALOG == "combined-v2":
    BASE = os.path.join(HERE, "results", "combined_usgs_isc_v2")
else:
    BASE = os.path.join(HERE, "results", "combined_usgs_isc")
OUT = os.path.join(BASE, CASE, "depth_clusters")
CASE_LABEL = (f"Catalog: {CATALOG} | M≥{MIN_MAG:g} input floor | "
              f"{'including' if INCLUDE_FIXED else 'excluding'} fixed depths | all years")
SOLVED_DEPTH = ("from location", "constrained by depth phases")
DEEP_MIN, DEEP_MAX = 40, 150
EPS_KM, MIN_EVENTS = 30.0, 4          # DBSCAN neighbourhood radius and minimum cluster size
MC_CUT = 4.4                          # reporting threshold; b is reported only for one type/source
BOX = dict(lat=(31.0, 35.0), lon=(74.0, 79.5))
GROUPS = [("intermediate", "#c62828"), ("uncertain", "#757575"), ("shallow", "#1565c0")]
os.makedirs(OUT, exist_ok=True)
matplotlib.rc("font", size=11)


def load():
    if CATALOG == "usgs":
        c = pd.read_csv(os.path.join(HERE, "data", "usgs_himachal_kashmir_M3plus_1970-2026.csv"))
        c = c.merge(pd.read_csv(os.path.join(HERE, "data", "usgs_depth_type.csv"))[["id", "depth_type"]],
                    on="id", how="left")
        c["source"] = "USGS"
        c["mag_source"] = "USGS"
    elif CATALOG == "combined-v2":
        c = pd.read_csv(os.path.join(HERE, "data", "combined_usgs_isc_catalog_v2.csv"))
    else:
        c = pd.read_csv(os.path.join(HERE, "data", "combined_usgs_isc_catalog.csv"))
    c["depth_type"] = c.depth_type.fillna("unknown")
    c["time"] = pd.to_datetime(c.time, utc=True, format="ISO8601").dt.tz_localize(None)
    c = c[c.latitude.between(*BOX["lat"]) & c.longitude.between(*BOX["lon"]) & c.mag.notna()]
    n_all = len(c)
    n_below_min = int((c.mag < MIN_MAG).sum())
    c = c[c.mag >= MIN_MAG].copy()
    if not INCLUDE_FIXED:
        c = c[c.depth_type.isin(SOLVED_DEPTH)]
    top, bottom = c.depth - c.depthError, c.depth + c.depthError
    c = c.assign(group=np.select([(top >= DEEP_MIN) & (c.depth <= DEEP_MAX), bottom < DEEP_MIN],
                                 ["intermediate", "shallow"], default="uncertain"))
    return c.sort_values("time").reset_index(drop=True), n_all, n_below_min


def b_mle(m, mc, dm=0.1):
    m = m[m >= mc]
    if len(m) < 5:
        return np.nan, np.nan, len(m)
    b = np.log10(np.e) / (m.mean() - (mc - dm / 2))
    return b, b / np.sqrt(len(m)), len(m)


def km_xy(lon, lat, lon0, lat0):
    return (lon - lon0) * 111.32 * np.cos(np.radians(lat0)), (lat - lat0) * 110.57


cat, n_all, n_below_min = load()
print(f"{CASE_LABEL}: {n_all} magnitude-bearing events in box, {n_below_min} below floor, {len(cat)} used after depth case, "
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
        scale_keys = set(zip(cl.magType.fillna("unknown"), cl.mag_source.fillna("unknown")))
        if len(scale_keys) == 1:
            bm, bs, nb = b_mle(cl.mag.values, MC_CUT)
            b_note = "single exact magType and mag_source"
        else:
            bm, bs, nb = np.nan, np.nan, int((cl.mag >= MC_CUT).sum())
            b_note = "not calculated: cluster mixes magnitude types or sources"
        rows.append(dict(depth_group=g, cluster=k, n_events=len(cl), n_M_ge_Mc=int((cl.mag >= MC_CUT).sum()),
                         lat_centre=round(lat0, 3), lon_centre=round(lon0, 3), max_extent_km=round(span, 1),
                         depth_median_km=round(cl.depth.median(), 1), depth_min_km=cl.depth.min(),
                         depth_max_km=cl.depth.max(), depthError_median_km=round(cl.depthError.median(), 1),
                         first=cl.time.min().date(), last=cl.time.max().date(),
                         mag_min=cl.mag.min(), mag_max=cl.mag.max(),
                         magTypes=";".join(sorted(cl.magType.fillna("unknown").astype(str).unique())),
                         magSources=";".join(sorted(cl.mag_source.fillna("unknown").astype(str).unique())),
                         b_mle=round(bm, 2) if bm == bm else None, b_err=round(bs, 2) if bs == bs else None,
                         b_mle_note=b_note,
                         interevent_CV=round(gaps.std() / gaps.mean(), 2) if len(gaps) > 2 else None,
                         share_fixed_or_unflagged=round((~cl.depth_type.isin(SOLVED_DEPTH)).mean(), 2),
                         sources=";".join(sorted(cl.source.astype(str).unique()))))
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, "depth_clusters.csv"), index=False)
cat.to_csv(os.path.join(OUT, "events_with_clusters.csv"), index=False)

# ----------------------------- MAP FIGURE ---------------------------------- #
fig, axs = plt.subplots(1, 3, figsize=(17, 6.2), sharex=True, sharey=True)
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
             "shallow": "Shallow: depth + error < 40 km"}[g]
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
az = np.radians(127.0)                                  # regional arc-parallel profile (see him_ides_himachal.py)
cat["along_km"] = x * np.sin(az) + y * np.cos(az)
for g, col in GROUPS:
    sub = cat[cat.group == g]
    inc = sub.cluster > 0
    known_error = sub.depthError.notna()
    show_error = inc & known_error
    missing_error_clustered = inc & ~known_error
    ax.errorbar(sub.along_km[show_error], sub.depth[show_error], yerr=sub.depthError[show_error], fmt="o", ms=5,
                color=col, ecolor=col, elinewidth=0.7, alpha=0.85, label=f"{g}, in a cluster")
    ax.scatter(sub.along_km[~inc], sub.depth[~inc], s=12, facecolors="none", edgecolors=col, lw=0.6, alpha=0.6,
               label=f"{g}, not clustered")
    if missing_error_clustered.any():
        ax.scatter(sub.along_km[missing_error_clustered], sub.depth[missing_error_clustered],
                   s=24, marker="^", color=col, edgecolors="k", linewidths=0.3,
                   label=f"{g}, clustered with no depth error")
ax.axhline(DEEP_MIN, ls=":", c="k", lw=0.8)
ax.set_ylim(DEEP_MAX, 0); ax.set_ylabel("Depth (km)"); ax.set_xlabel("Distance along the arc, azimuth 127° (km)")
ax.legend(fontsize=8, ncol=3, loc="lower left"); ax.grid(alpha=0.2)
clipped_errors = int(((cat.depth - cat.depthError < 0) | (cat.depth + cat.depthError > DEEP_MAX)).fillna(False).sum())
ax.set_title(f"Events projected onto one arc-parallel line; error bars are source-reported   ·   {CASE_LABEL}", fontsize=10)
fig.text(0.01, 0.01, f"Depth errors are clipped at 0–{DEEP_MAX} km (N extending outside={clipped_errors}); missing errors are marked separately.",
         fontsize=8, color="0.35")
plt.tight_layout(); fig.savefig(os.path.join(OUT, "DepthClusters_section.png"), dpi=200); plt.close(fig)

pd.set_option("display.width", 250)
print(tab.drop(columns=["sources"]).to_string(index=False) if len(tab) else "no clusters")
print("Outputs in", OUT)
