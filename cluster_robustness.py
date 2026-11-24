"""Robustness of the intermediate-depth clusters: DBSCAN settings grid and Monte Carlo over location errors."""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from sklearn.cluster import DBSCAN

from seis_utils import (BOX, CLUSTER_POLY, DEEP_MAX, DEEP_MIN, EARTH_R, HERE, load_catalog)

OUT = os.path.join(HERE, "results", "robustness")
EPS_GRID, MIN_GRID = [15, 20, 30, 45], [3, 4, 5, 6]
EPS_KM, MIN_EVENTS = 30.0, 4
N_MC, SEED = 500, 0
os.makedirs(OUT, exist_ok=True)


def dbscan(lat, lon, eps_km, min_events):
    """DBSCAN labels on a haversine metric; all -1 if too few events."""
    if len(lat) < min_events:
        return np.full(len(lat), -1)
    X = np.radians(np.c_[lat, lon])
    return DBSCAN(eps=eps_km / EARTH_R, min_samples=min_events, metric="haversine").fit_predict(X)


def band_share(labels, inside):
    """Largest share of the in-outline events that sit in one cluster."""
    if inside.sum() == 0:
        return np.nan
    lab = labels[inside]
    lab = lab[lab >= 0]
    return 0.0 if len(lab) == 0 else np.bincount(lab).max() / inside.sum()


def grid(cat_name, df):
    """Cluster the intermediate events for every radius / minimum-size setting."""
    deep = df[df.depth_class == "intermediate"]
    inside = Path(CLUSTER_POLY).contains_points(deep[["longitude", "latitude"]].values)
    rows = []
    for eps in EPS_GRID:
        for mn in MIN_GRID:
            lab = dbscan(deep.latitude.values, deep.longitude.values, eps, mn)
            sizes = np.bincount(lab[lab >= 0]) if (lab >= 0).any() else np.array([], int)
            share = band_share(lab, inside)
            rows.append(dict(catalog=cat_name, eps_km=eps, min_events=mn, n_intermediate=len(deep),
                             n_in_outline=int(inside.sum()), n_clusters=len(sizes),
                             n_clustered=int(sizes.sum()), largest_cluster=int(sizes.max()) if len(sizes) else 0,
                             band_share_in_one_cluster=share, band_recovered=bool(share >= 0.5)))
    return pd.DataFrame(rows)


def monte_carlo(cat_name, df, horiz_default):
    """Resample depths and epicentres within their errors N_MC times.
    Returns per-event probabilities and per-run cluster results."""
    rng = np.random.default_rng(SEED)
    z, ez = df.depth.values, df.depthError.fillna(0).values
    has_h = df.horizontalError.notna().values
    eh = np.where(has_h, df.horizontalError.values, horiz_default)
    n = len(df)
    deep = (df.depth_class == "intermediate").values
    lat, lon = df.latitude.values[deep], df.longitude.values[deep]
    inside = Path(CLUSTER_POLY).contains_points(np.c_[lon, lat])
    p_int, p_clu = np.zeros(n), np.zeros(n)
    runs = []
    for k in range(N_MC):
        zz = z + rng.normal(0, 1, n) * ez
        p_int += (zz >= DEEP_MIN) & (zz <= DEEP_MAX)
        dx, dy = rng.normal(0, 1, deep.sum()) * eh[deep], rng.normal(0, 1, deep.sum()) * eh[deep]
        la = lat + dy / 110.57
        lo = lon + dx / (111.32 * np.cos(np.radians(lat)))
        lab = dbscan(la, lo, EPS_KM, MIN_EVENTS)
        p_clu[deep] += lab >= 0
        share = band_share(lab, inside)
        runs.append(dict(catalog=cat_name, run=k, n_intermediate=int(deep.sum()),
                         n_clusters=int(lab.max() + 1) if (lab >= 0).any() else 0,
                         n_in_outline=int(inside.sum()), band_share_in_one_cluster=share,
                         band_recovered=bool(share >= 0.5)))
    ev = df[["id", "time", "latitude", "longitude", "depth", "depthError", "depth_class", "mag", "magType"]].copy()
    ev["horizontal_error_used_km"] = eh
    ev["horizontal_error_assumed"] = ~has_h
    ev["p_intermediate"] = p_int / N_MC
    ev["p_in_cluster"] = np.where(deep, p_clu / N_MC, np.nan)
    ev.to_csv(os.path.join(OUT, f"mc_event_probabilities_{cat_name}.csv"), index=False)
    return ev, pd.DataFrame(runs)


def main():
    """Run the settings grid and the Monte Carlo for both catalogs, write tables and figures."""
    usgs, _ = load_catalog("usgs", include_fixed=False, verbose=False)
    comb, _ = load_catalog("combined", include_fixed=False, verbose=False)
    horiz_default = float(usgs.horizontalError.median())
    print(f"Horizontal error assumed where missing: {horiz_default:.1f} km "
          f"(median USGS horizontalError, solved depths)")
    g = pd.concat([grid("usgs", usgs), grid("combined", comb)], ignore_index=True)
    g.to_csv(os.path.join(OUT, "dbscan_grid.csv"), index=False)

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, cat in zip(axs, ("usgs", "combined")):
        d = g[g.catalog == cat].pivot(index="min_events", columns="eps_km", values="band_share_in_one_cluster")
        nc = g[g.catalog == cat].pivot(index="min_events", columns="eps_km", values="n_clusters")
        im = ax.imshow(d.values, cmap="Blues", vmin=0, vmax=1, origin="lower", aspect="auto")
        for i in range(d.shape[0]):
            for j in range(d.shape[1]):
                v = d.values[i, j]
                ax.text(j, i, f"{v:.0%}\n{int(nc.values[i, j])} cl.", ha="center", va="center", fontsize=8,
                        color="white" if v > 0.6 else "black")
        ax.set_xticks(range(d.shape[1]), d.columns); ax.set_yticks(range(d.shape[0]), d.index)
        ax.set_xlabel("DBSCAN radius (km)"); ax.set_ylabel("Minimum events per cluster")
        n_in = int(g[g.catalog == cat].n_in_outline.iloc[0])
        ax.set_title(f"{'USGS' if cat == 'usgs' else 'USGS + ISC'} ({n_in} events in the outline)", fontsize=10)
    fig.colorbar(im, ax=axs, shrink=0.8, label="Share in one cluster")
    fig.suptitle("Share of the Kishtwar-outline intermediate events that DBSCAN puts in one cluster "
                 "(cl. = clusters found in the whole region)", fontsize=10, y=1.02)
    fig.savefig(os.path.join(OUT, "dbscan_grid.png"), dpi=200, bbox_inches="tight"); plt.close(fig)

    evs, runs = {}, []
    for cat, df in (("usgs", usgs), ("combined", comb)):
        evs[cat], r = monte_carlo(cat, df, horiz_default)
        runs.append(r)
    runs = pd.concat(runs, ignore_index=True)
    summ = runs.groupby("catalog").agg(runs=("run", "size"), n_intermediate=("n_intermediate", "first"),
                                       n_in_outline=("n_in_outline", "first"),
                                       clusters_median=("n_clusters", "median"),
                                       clusters_p5=("n_clusters", lambda x: np.percentile(x, 5)),
                                       clusters_p95=("n_clusters", lambda x: np.percentile(x, 95)),
                                       band_recovered_share=("band_recovered", "mean"),
                                       band_share_median=("band_share_in_one_cluster", "median")).reset_index()
    summ.to_csv(os.path.join(OUT, "mc_runs_summary.csv"), index=False)

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    poly = np.array(CLUSTER_POLY + [CLUSTER_POLY[0]])
    for ax, cat in zip(axs, ("usgs", "combined")):
        e = evs[cat]
        show = e[e.p_intermediate > 0].sort_values("p_intermediate")
        deep = e[e.depth_class == "intermediate"]
        ax.scatter(e.longitude, e.latitude, s=3, c="0.85", lw=0)
        sc = ax.scatter(show.longitude, show.latitude, c=show.p_intermediate, cmap="Reds", vmin=0, vmax=1,
                        s=10, edgecolors="none")
        ax.scatter(deep.longitude, deep.latitude, s=15 + 80 * deep.p_in_cluster, facecolors="none",
                   edgecolors="k", linewidths=0.8)
        ax.plot(poly[:, 0], poly[:, 1], "k-", lw=0.8)
        ax.set_xlim(*BOX["lon"]); ax.set_ylim(*BOX["lat"]); ax.set_aspect(1 / np.cos(np.radians(33)))
        ax.set_xlabel("Longitude (°)"); ax.grid(alpha=0.2)
        s = summ[summ.catalog == cat].iloc[0]
        ax.set_title(f"{'USGS' if cat == 'usgs' else 'USGS + ISC'}, solved depths, {N_MC} resamplings\n"
                     f"Kishtwar band ({int(s.n_in_outline)} events) in one cluster in "
                     f"{s.band_recovered_share:.0%} of runs", fontsize=10)
    axs[0].set_ylabel("Latitude (°)")
    fig.colorbar(sc, ax=axs, shrink=0.8, label="Probability of intermediate depth (40-150 km)")
    fig.text(0.01, 0.01, "Fill: probability of 40-150 km from depth errors. Black rings: rule-based intermediate "
             "events, ring size = probability of staying in a DBSCAN cluster (30 km, 4 events) when epicentres "
             f"are shifted within their errors. Missing horizontal errors set to {horiz_default:.1f} km (assumption).",
             fontsize=8)
    fig.savefig(os.path.join(OUT, "mc_probability_map.png"), dpi=200, bbox_inches="tight"); plt.close(fig)

    pd.set_option("display.width", 200)
    print(g.round(2).to_string(index=False))
    print(summ.round(2).to_string(index=False))
    for cat in evs:
        e = evs[cat]
        d = e[e.depth_class == "intermediate"]
        print(f"{cat}: events with p_intermediate >= 0.9: {(e.p_intermediate >= 0.9).sum()}, "
              f">= 0.5: {(e.p_intermediate >= 0.5).sum()}; rule-based intermediate: {len(d)} "
              f"(min p_intermediate among them {d.p_intermediate.min():.2f}); "
              f"horizontal error assumed for {int(e.horizontal_error_assumed.sum())} of {len(e)} events")
    print("Outputs in", OUT)


if __name__ == "__main__":
    main()
