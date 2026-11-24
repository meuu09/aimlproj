"""Region-wide candidate-sequence and migration search with window sensitivity (results/migration/)."""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path

from seis_utils import (CLUSTER_POLY, HERE, SOLVED_DEPTH, candidate_sequences, hypo_dist_km, load_catalog,
                        migration_fit, resolve_overlaps)

OUT = os.path.join(HERE, "results", "migration")
CATALOGS = [("usgs", "usgs", False), ("combined", "combined", False), ("combined_nodup", "combined", True)]
DAYS, KMS = [5, 10, 20, 30, 60], [10, 20, 30, 50]
PAPER = (20, 20)
MIN_EVENTS = 3
os.makedirs(OUT, exist_ok=True)


def populations(df):
    """Event sets searched: Kishtwar outline (intermediate), all intermediate, all solved depths."""
    inside = Path(CLUSTER_POLY).contains_points(df[["longitude", "latitude"]].values)
    return {"Kishtwar outline, intermediate": df[inside & (df.depth_class == "intermediate")],
            "region, intermediate": df[df.depth_class == "intermediate"],
            "region, all depth classes": df}


def describe(d, c):
    """Fits and quality numbers for one candidate."""
    g = d.loc[c["members"]].sort_values("time")
    e = d.loc[c["first"]]
    t = (g.time - e.time).dt.total_seconds().values / 86400
    r = hypo_dist_km(e.latitude, e.longitude, e.depth, g.latitude.values, g.longitude.values, g.depth.values)
    mm = np.sort(g.mag.values)[::-1]
    err = np.sqrt(np.nanmedian(g.horizontalError) ** 2 + np.nanmedian(g.depthError) ** 2)
    row = dict(start=e.time, n_events=len(g), duration_days=t.max(), max_dist_km=r.max(),
               depth_classes=";".join(g.depth_class), mags=", ".join(f"{m:.1f}" for m in g.mag),
               dM_largest_minus_second=mm[0] - mm[1], first_is_largest=bool(g.mag.iloc[0] == mm[0]),
               median_location_error_km=err, resolvable=bool(r.max() > err) if np.isfinite(err) else False,
               ids=";".join(g.id))
    row.update(migration_fit(t, r))
    return row, t, r


def main():
    """Grid search, candidate table, figure of candidate sequences."""
    grid, cands_out, plots = [], [], []
    for name, src, nodup in CATALOGS:
        for case in ("excluding_fixed", "including_fixed"):
            df, _ = load_catalog(src, True, None, nodup, verbose=False)
            if case == "excluding_fixed":
                df = df[df.depth_type.isin(SOLVED_DEPTH)]
            for pop, sub in populations(df).items():
                for days in DAYS:
                    for km in KMS:
                        d, cands = candidate_sequences(sub, days, km, MIN_EVENTS)
                        kept = resolve_overlaps(cands)
                        grid.append(dict(catalog=name, case=case, population=pop, days=days, km=km,
                                         n_events=len(sub), n_candidates=len(cands), n_after_overlap_rule=len(kept)))
                        if (days, km) != PAPER:
                            continue
                        keep_ids = {c["ids"] for c in kept}
                        for c in cands:
                            row, t, r = describe(d, c)
                            row.update(catalog=name, case=case, population=pop, days=days, km=km,
                                       kept_after_overlap_rule=c["ids"] in keep_ids)
                            cands_out.append(row)
                            if c["ids"] in keep_ids and pop != "region, all depth classes" and name == "combined":
                                plots.append((f"{name} {case} {row['start']:%Y-%m-%d}", t, r, row))
    g = pd.DataFrame(grid)
    g.to_csv(os.path.join(OUT, "window_sensitivity.csv"), index=False)
    c = pd.DataFrame(cands_out)
    c.to_csv(os.path.join(OUT, "sequence_candidates_20d_20km.csv"), index=False)
    if plots:
        fig, axs = plt.subplots(1, len(plots), figsize=(4 * len(plots), 3.6), squeeze=False)
        for ax, (title, t, r, row) in zip(axs[0], plots):
            ax.scatter(t, r, c="k", s=25)
            tt = np.linspace(0, max(t.max(), 1e-3), 100)
            if "D_fit_m2_s" in row:
                ax.plot(tt, np.sqrt(4 * np.pi * row["D_fit_m2_s"] * tt * 86400) / 1000, "b--", label="diffusive fit")
            if "v_origin_km_day" in row:
                ax.plot(tt, row["v_origin_km_day"] * tt, "r-", label="linear through origin")
            ax.axhline(row["median_location_error_km"], color="0.5", ls=":", label="median location error")
            ax.set_title(title, fontsize=8); ax.set_xlabel("Days after first event"); ax.grid(alpha=0.2)
        axs[0][0].set_ylabel("3-D distance from first event (km)"); axs[0][0].legend(fontsize=7)
        plt.tight_layout(); fig.savefig(os.path.join(OUT, "sequence_candidates.png"), dpi=200); plt.close(fig)
    pd.set_option("display.width", 220)
    p = g[(g.days == PAPER[0]) & (g.km == PAPER[1])]
    print(p.to_string(index=False))
    if len(c):
        cols = ["catalog", "case", "population", "start", "n_events", "max_dist_km", "median_location_error_km",
                "resolvable", "dM_largest_minus_second", "v_free_km_day", "trend_perm_p", "kept_after_overlap_rule"]
        print(c[[k for k in cols if k in c]].round(2).to_string(index=False))
    print("Outputs in", OUT)


if __name__ == "__main__":
    main()
