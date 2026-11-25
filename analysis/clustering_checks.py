"""3-D vs 2-D DBSCAN, k-distance, polygon bias, PCA strikes, space-time feasibility and base map (results/clustering/)."""
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from shapely.geometry import Polygon
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

from seis_utils import BOX, CLUSTER_POLY, DATA, EARTH_R, HERE, load_catalog, xy_km

OUT = os.path.join(HERE, "results", "clustering")
EPS_KM, MIN_EVENTS = 30.0, 4
LAT0, LON0 = np.mean(BOX["lat"]), np.mean(BOX["lon"])
TOWNS = {"Srinagar": (74.80, 34.08), "Jammu": (74.86, 32.73), "Kishtwar": (75.77, 33.31), "Chamba": (76.13, 32.55),
         "Dharamshala": (76.32, 32.22), "Shimla": (77.17, 31.10), "Leh": (77.58, 34.16)}
os.makedirs(OUT, exist_ok=True)


def labels_2d(d, eps=EPS_KM, mn=MIN_EVENTS):
    """Haversine DBSCAN on [lat, lon] in radians."""
    if len(d) < mn:
        return np.full(len(d), -1)
    return DBSCAN(eps=eps / EARTH_R, min_samples=mn, metric="haversine").fit_predict(
        np.radians(d[["latitude", "longitude"]].values))


def coords_3d(d):
    """East, north (local tangent plane around the box centre) and depth, all in km."""
    x, y = xy_km(d.longitude.values, d.latitude.values, LON0, LAT0)
    return np.c_[x, y, d.depth.values]


def labels_3d(d, eps=EPS_KM, mn=MIN_EVENTS):
    """Euclidean DBSCAN on km coordinates including depth."""
    if len(d) < mn:
        return np.full(len(d), -1)
    return DBSCAN(eps=eps, min_samples=mn).fit_predict(coords_3d(d))


def pca_azimuth(d):
    """Azimuth (0-180 deg, clockwise from north) of the first principal axis of epicentres in km."""
    if len(d) < 3:
        return np.nan, np.nan
    x, y = xy_km(d.longitude.values, d.latitude.values, LON0, LAT0)
    xy = np.c_[x - x.mean(), y - y.mean()]
    _, s, vt = np.linalg.svd(xy, full_matrices=False)
    return float(np.degrees(np.arctan2(vt[0][0], vt[0][1])) % 180), float(s[0] ** 2 / (s ** 2).sum())


def compare(name, d):
    """2-D vs 3-D clusters and the Kishtwar band, for one catalog's intermediate events."""
    l2, l3 = labels_2d(d), labels_3d(d)
    inside = Path(CLUSTER_POLY).contains_points(d[["longitude", "latitude"]].values)
    rows = []
    for dim, lab in (("2-D", l2), ("3-D", l3)):
        for k in sorted(set(lab[lab >= 0])):
            m = lab == k
            az, var = pca_azimuth(d[m])
            rows.append(dict(catalog=name, dbscan=dim, cluster=int(k), n=int(m.sum()),
                             lat=d.latitude[m].mean(), lon=d.longitude[m].mean(), depth_median=d.depth[m].median(),
                             depth_min=d.depth[m].min(), depth_max=d.depth[m].max(),
                             n_in_kishtwar_outline=int((m & inside).sum()), pca_azimuth=az, pca_variance_share=var))
    summary = dict(catalog=name, n_intermediate=len(d), n_clusters_2d=len(set(l2[l2 >= 0])),
                   n_clusters_3d=len(set(l3[l3 >= 0])), noise_2d=float((l2 < 0).mean()),
                   noise_3d=float((l3 < 0).mean()), adjusted_rand_2d_vs_3d=float(adjusted_rand_score(l2, l3)),
                   n_in_outline=int(inside.sum()))
    for dim, lab in (("2d", l2), ("3d", l3)):
        b = lab[inside & (lab >= 0)]
        best = np.bincount(b).argmax() if len(b) else -1
        clu = lab == best if best >= 0 else np.zeros(len(d), bool)
        summary[f"band_share_in_one_cluster_{dim}"] = (clu & inside).sum() / max(inside.sum(), 1)
        summary[f"jaccard_outline_vs_cluster_{dim}"] = (clu & inside).sum() / max((clu | inside).sum(), 1)
    return pd.DataFrame(rows), summary, l2, l3


def polygon_variants(name, d, l2):
    """Membership of the hand-drawn outline vs buffered / shrunk outlines and the DBSCAN band cluster."""
    base = Polygon(CLUSTER_POLY)
    rows = []
    for buf in (-0.2, -0.1, 0.0, 0.1, 0.2):
        poly = base.buffer(buf, join_style=2)
        inside = Path(np.array(poly.exterior.coords)).contains_points(d[["longitude", "latitude"]].values)
        b = l2[inside & (l2 >= 0)]
        best = np.bincount(b).argmax() if len(b) else -1
        rows.append(dict(catalog=name, buffer_deg=buf, n_inside=int(inside.sum()),
                         n_in_largest_cluster=int((l2[inside] == best).sum()) if best >= 0 else 0,
                         depth_median=d.depth[inside].median() if inside.any() else np.nan))
    return rows


def st_feasibility(name, d, km=30.0, days=365.0):
    """How many intermediate events have >= MIN_EVENTS-1 neighbours within km and days (ST-DBSCAN core points)."""
    xyz = coords_3d(d)[:, :2]
    t = d.time.values.astype("datetime64[s]").astype(float) / 86400
    near = (np.hypot(*(xyz[:, None, :] - xyz[None, :, :]).transpose(2, 0, 1)) <= km) & (
        np.abs(t[:, None] - t[None, :]) <= days)
    cores = int(((near.sum(1) - 1) >= MIN_EVENTS - 1).sum())
    return dict(catalog=name, n_intermediate=len(d), space_km=km, time_days=days, n_core_points=cores)


def kdist_plot(sets):
    """Sorted distance to the 3rd nearest neighbour (k = MIN_EVENTS - 1), 2-D and 3-D."""
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for name, d in sets:
        for ax, X, metric in ((axs[0], np.radians(d[["latitude", "longitude"]].values), "haversine"),
                              (axs[1], coords_3d(d), "euclidean")):
            nn = NearestNeighbors(n_neighbors=MIN_EVENTS, metric=metric).fit(X)
            dist = nn.kneighbors(X)[0][:, -1] * (EARTH_R if metric == "haversine" else 1)
            ax.plot(np.sort(dist)[::-1], label=name)
    for ax, t in zip(axs, ("2-D epicentral", "3-D hypocentral")):
        ax.axhline(EPS_KM, ls=":", c="k", lw=0.8)
        ax.set_title(f"{t}: distance to {MIN_EVENTS - 1}rd neighbour", fontsize=10)
        ax.set_xlabel("Events, sorted"); ax.set_ylabel("km"); ax.grid(alpha=0.2); ax.set_yscale("log")
    axs[0].legend(fontsize=8)
    plt.tight_layout(); fig.savefig(os.path.join(OUT, "k_distance.png"), dpi=200); plt.close(fig)


def gaf534_azimuth():
    """Trace azimuth of GEM fault GAF_534 (reverse, Himalayan front in Himachal), an independent strike."""
    g = json.load(open(os.path.join(DATA, "gem_active_faults_box.geojson")))
    for f in g["features"]:
        if f["properties"].get("catalog_id") == "GAF_534":
            c = np.array(f["geometry"]["coordinates"] if f["geometry"]["type"] == "LineString"
                         else sum(f["geometry"]["coordinates"], []))
            return pca_azimuth(pd.DataFrame({"longitude": c[:, 0], "latitude": c[:, 1]}))[0]
    return np.nan


def base_map(sets, fm):
    """Epicentres by depth class, GEM faults, Kishtwar outline, GCMT mechanisms and reference towns."""
    g = json.load(open(os.path.join(DATA, "gem_active_faults_box.geojson")))
    fig, axs = plt.subplots(1, len(sets), figsize=(8 * len(sets), 6.8), squeeze=False)
    col = {"thrust": "#d62728", "thrust with strike-slip": "#ff7f0e", "strike-slip": "#2ca02c",
           "normal": "#1f77b4", "normal with strike-slip": "#17becf", "unknown (oblique)": "#7f7f7f"}
    for ax, (name, d) in zip(axs[0], sets):
        for f in g["features"]:
            lines = [f["geometry"]["coordinates"]] if f["geometry"]["type"] == "LineString" else f["geometry"]["coordinates"]
            for l in lines:
                l = np.array(l)
                ax.plot(l[:, 0], l[:, 1], c="0.35", lw=1)
        for cls, c, s in (("shallow", "#9ecae1", 6), ("uncertain", "#bdbdbd", 6), ("intermediate", "#c62828", 22),
                          ("deep", "k", 22)):
            e = d[d.depth_class == cls]
            ax.scatter(e.longitude, e.latitude, s=s, c=c, lw=0, label=f"{cls} ({len(e)})", zorder=3)
        p = np.array(CLUSTER_POLY + [CLUSTER_POLY[0]])
        ax.plot(p[:, 0], p[:, 1], "k--", lw=1, label="Kishtwar outline")
        for reg, gg in fm.groupby("regime"):
            ax.scatter(gg.centroid_longitude, gg.centroid_latitude, marker="*", s=140, c=col.get(reg, "k"),
                       edgecolors="k", lw=0.5, zorder=4, label=f"GCMT {reg} ({len(gg)})")
        for town, (lo, la) in TOWNS.items():
            ax.plot(lo, la, "ks", ms=3)
            ax.annotate(town, (lo, la), xytext=(3, 3), textcoords="offset points", fontsize=8)
        ax.set_xlim(*BOX["lon"]); ax.set_ylim(*BOX["lat"]); ax.set_aspect(1 / np.cos(np.radians(LAT0)))
        ax.set_title(f"{name}, solved depths (all years)", fontsize=10)
        ax.set_xlabel("Longitude (°)"); ax.set_ylabel("Latitude (°)"); ax.legend(fontsize=7, loc="lower left")
    fig.text(0.01, 0.005, "Grey lines: GEM Global Active Faults (Styron & Pagani 2020, CC BY-SA 4.0); coverage of the "
             "Indian NW Himalaya is incomplete (no frontal-thrust traces near Kishtwar). Stars: Global CMT centroids, "
             "regime from axis plunges (Zoback 1992). Town positions approximate.", fontsize=8)
    fig.savefig(os.path.join(OUT, "base_map.png"), dpi=200, bbox_inches="tight"); plt.close(fig)


def main():
    """Run every check for the USGS and combined catalogs (solved depths, all years)."""
    sets = []
    for name, src, nodup in (("usgs", "usgs", False), ("combined", "combined", False),
                             ("combined_nodup", "combined", True)):
        df, _ = load_catalog(src, False, None, nodup, verbose=False)
        sets.append((name, df))
    clusters, summ, poly, st, pca = [], [], [], [], []
    for name, df in sets:
        d = df[df.depth_class == "intermediate"].reset_index(drop=True)
        c, s, l2, l3 = compare(name, d)
        clusters.append(c); summ.append(s)
        poly += polygon_variants(name, d, l2)
        st += [st_feasibility(name, d, km, days) for km in (20, 30, 50) for days in (30, 365, 1825)]
        inside = Path(CLUSTER_POLY).contains_points(df[["longitude", "latitude"]].values)
        for lab, sub in (("all events M>=4.4", df[df.mag >= 4.4]), ("intermediate", d),
                         ("Kishtwar outline, all", df[inside]),
                         ("largest 2-D intermediate cluster", d[l2 == np.bincount(l2[l2 >= 0]).argmax()]
                          if (l2 >= 0).any() else d.iloc[:0])):
            az, var = pca_azimuth(sub)
            pca.append(dict(catalog=name, events=lab, n=len(sub), pca_azimuth=az, variance_share=var))
    pca.append(dict(catalog="GEM fault GAF_534 trace (Himalayan front, Himachal)", events="fault trace", n=np.nan,
                    pca_azimuth=gaf534_azimuth(), variance_share=np.nan))
    pd.concat(clusters).to_csv(os.path.join(OUT, "dbscan_2d_3d_clusters.csv"), index=False)
    pd.DataFrame(summ).to_csv(os.path.join(OUT, "dbscan_2d_vs_3d_summary.csv"), index=False)
    pd.DataFrame(poly).to_csv(os.path.join(OUT, "polygon_variants.csv"), index=False)
    pd.DataFrame(st).to_csv(os.path.join(OUT, "space_time_feasibility.csv"), index=False)
    pd.DataFrame(pca).to_csv(os.path.join(OUT, "pca_strikes.csv"), index=False)
    kdist_plot([(n, d[d.depth_class == "intermediate"]) for n, d in sets[:2]])
    fm_path = os.path.join(HERE, "results", "geology", "focal_mechanisms.csv")
    fm = pd.read_csv(fm_path) if os.path.exists(fm_path) else pd.DataFrame(columns=["regime"])
    base_map(sets[:2], fm)
    pd.set_option("display.width", 220)
    for t in (pd.DataFrame(summ), pd.concat(clusters), pd.DataFrame(poly), pd.DataFrame(pca), pd.DataFrame(st)):
        print(t.round(2).to_string(index=False))
    print("Outputs in", OUT)


if __name__ == "__main__":
    main()
