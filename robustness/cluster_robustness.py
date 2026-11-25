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


SCENARIOS = [("errors as 1 sigma", 1.0, 1.0), ("errors as 90 % bounds", 1 / 1.645, 1.0),
             ("errors understated x2", 2.0, 1.0), ("missing horizontal error = 2 x median", 1.0, 2.0)]
THRESHOLDS = (0.5, 0.7, 0.8, 0.9)


def monte_carlo(cat_name, df, horiz_default, sigma_scale=1.0, h_fill_scale=1.0, n_runs=N_MC, seed=SEED):
    """Joint resampling: perturb depth (only where an error exists) and epicentre, reclassify each event,
    rerun DBSCAN on the realised intermediate events and score the Kishtwar band against its catalog members."""
    rng = np.random.default_rng(seed)
    n = len(df)
    z, ez = df.depth.values, df.depthError.values
    has_z = np.isfinite(ez)
    has_h = df.horizontalError.notna().values
    eh = np.where(has_h, df.horizontalError.values, horiz_default * h_fill_scale) * sigma_scale
    lat0, lon0 = df.latitude.values, df.longitude.values
    orig_int = (df.depth_class == "intermediate").values
    band = orig_int & Path(CLUSTER_POLY).contains_points(np.c_[lon0, lat0])
    p_int, p_clu = np.zeros(n), np.zeros(n)
    runs = []
    for k in range(n_runs):
        zz = np.where(has_z, z + rng.normal(0, 1, n) * np.nan_to_num(ez) * sigma_scale, z)
        real = has_z & (zz >= DEEP_MIN) & (zz <= DEEP_MAX)
        p_int += real
        la = lat0 + rng.normal(0, 1, n) * eh / 110.57
        lo = lon0 + rng.normal(0, 1, n) * eh / (111.32 * np.cos(np.radians(lat0)))
        lab = np.full(n, -1)
        lab[real] = dbscan(la[real], lo[real], EPS_KM, MIN_EVENTS)
        p_clu += lab >= 0
        bl = lab[band]
        bl = bl[bl >= 0]
        share = np.bincount(bl).max() / band.sum() if len(bl) else 0.0
        best = np.bincount(bl).argmax() if len(bl) else -1
        size = int((lab == best).sum()) if best >= 0 else 0
        r = dict(catalog=cat_name, run=k, n_realised_intermediate=int(real.sum()), n_band=int(band.sum()),
                 n_clusters=int(len(np.unique(lab[lab >= 0]))), band_share_in_one_cluster=share,
                 contamination=(size - int((bl == best).sum())) / size if size else np.nan,
                 fragments=int(len(np.unique(bl))), band_dropped_from_intermediate=int((band & ~real).sum()))
        r.update({f"recovered_{int(t * 100)}": bool(share >= t) for t in THRESHOLDS})
        runs.append(r)
    ev = df[["id", "time", "latitude", "longitude", "depth", "depthError", "depth_class", "mag", "magType"]].copy()
    ev["horizontal_error_used_km"] = eh
    ev["horizontal_error_assumed"] = ~has_h
    ev["depth_error_missing"] = ~has_z
    ev["p_intermediate"] = np.where(has_z, p_int / n_runs, np.nan)
    ev["p_in_cluster"] = np.where(orig_int, p_clu / n_runs, np.nan)
    ev["in_kishtwar_band"] = band
    return ev, pd.DataFrame(runs)


def summarise_runs(runs):
    """Per catalog (and scenario): recovery at each threshold, contamination, fragmentation."""
    keys = [k for k in ("scenario", "catalog") if k in runs]
    agg = dict(runs=("run", "size"), n_band=("n_band", "first"),
               realised_intermediate_median=("n_realised_intermediate", "median"),
               clusters_median=("n_clusters", "median"), clusters_p5=("n_clusters", lambda x: np.percentile(x, 5)),
               clusters_p95=("n_clusters", lambda x: np.percentile(x, 95)),
               band_share_median=("band_share_in_one_cluster", "median"),
               contamination_median=("contamination", "median"), fragments_median=("fragments", "median"),
               band_dropped_median=("band_dropped_from_intermediate", "median"))
    agg.update({f"recovered_{int(t * 100)}_share": (f"recovered_{int(t * 100)}", "mean") for t in THRESHOLDS})
    return runs.groupby(keys).agg(**agg).reset_index()


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

    nodup, _ = load_catalog("combined", include_fixed=False, drop_duplicates=True, verbose=False)
    cats = (("usgs", usgs), ("combined", comb), ("combined_nodup", nodup))
    evs, runs, scen = {}, [], []
    for cat, df in cats:
        evs[cat], r = monte_carlo(cat, df, horiz_default)
        evs[cat].to_csv(os.path.join(OUT, f"mc_event_probabilities_{cat}.csv"), index=False)
        runs.append(r)
        for name, s, h in SCENARIOS:
            scen.append(monte_carlo(cat, df, horiz_default, s, h)[1].assign(scenario=name))
    runs = pd.concat(runs, ignore_index=True)
    summ = summarise_runs(runs)
    summ["band_recovered_share"] = summ.recovered_50_share
    summ.to_csv(os.path.join(OUT, "mc_runs_summary.csv"), index=False)
    ssum = summarise_runs(pd.concat(scen, ignore_index=True))
    ssum.to_csv(os.path.join(OUT, "mc_scenarios_summary.csv"), index=False)

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), sharey=True)
    poly = np.array(CLUSTER_POLY + [CLUSTER_POLY[0]])
    for ax, cat in zip(axs, ("usgs", "combined")):
        e = evs[cat]
        e = e.assign(p_intermediate=e.p_intermediate.fillna(0))
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
                     f"Kishtwar band ({int(s.n_band)} events): >=50 % in one cluster in "
                     f"{s.recovered_50_share:.0%} of runs, >=80 % in {s.recovered_80_share:.0%}", fontsize=10)
    axs[0].set_ylabel("Latitude (°)")
    fig.colorbar(sc, ax=axs, shrink=0.8, label="Probability of intermediate depth (40-150 km)")
    fig.text(0.01, 0.01, "Fill: share of realisations at 40-150 km (events without a depth error not shown). "
             "Black rings: rule-based intermediate events, ring size = share of joint realisations (depth and "
             "epicentre) in which they are intermediate and in a DBSCAN cluster (30 km, 4 events). Missing horizontal "
             f"errors set to {horiz_default:.1f} km (assumption). Percentages are simulation recovery rates, "
             "not probabilities that a structure exists.",
             fontsize=8)
    fig.savefig(os.path.join(OUT, "mc_probability_map.png"), dpi=200, bbox_inches="tight"); plt.close(fig)

    pd.set_option("display.width", 200)
    print(g.round(2).to_string(index=False))
    print(summ.round(2).to_string(index=False))
    print(ssum.round(2).to_string(index=False))
    for cat in evs:
        e = evs[cat]
        d = e[e.depth_class == "intermediate"]
        print(f"{cat}: events with p_intermediate >= 0.9: {(e.p_intermediate >= 0.9).sum()}, "
              f">= 0.5: {(e.p_intermediate >= 0.5).sum()}; rule-based intermediate: {len(d)} "
              f"(min p_intermediate among them {d.p_intermediate.min():.2f}); "
              f"no depth error: {int(e.depth_error_missing.sum())}; "
              f"horizontal error assumed for {int(e.horizontal_error_assumed.sum())} of {len(e)} events")
    print("Outputs in", OUT)


if __name__ == "__main__":
    main()
