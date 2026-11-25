"""Collect key numbers from results/ into one long CSV, for before/after comparison."""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(HERE, "results")
CATS = {"usgs": R, "combined": os.path.join(R, "combined_usgs_isc"), "combined_nodup": os.path.join(R, "combined_usgs_isc_nodup")}


def read(path):
    """CSV as DataFrame, or None if missing."""
    return pd.read_csv(path) if os.path.exists(path) else None


def main(out):
    """Write metric, catalog, case, value rows to `out`."""
    rows = []
    add = lambda m, c, k, v: rows.append(dict(metric=m, catalog=c, case=k, value=v))
    co = read(os.path.join(HERE, "data", "combined_usgs_isc_catalog.csv"))
    add("combined_rows", "combined", "", len(co))
    add("matched_usgs_isc", "combined", "", int((co.source == "USGS+ISC").sum()))
    add("possible_duplicates", "combined", "", int(co.possible_duplicate_of.notna().sum()))
    add("solved_depths", "combined", "", int(co.depth_type.isin(["from location", "constrained by depth phases"]).sum()))
    add("matches_over_50km", "combined", "", int((co.match_dist_km > 50).sum()))
    fun = read(os.path.join(R, "qc", "filter_funnel.csv"))
    for _, r in fun.iterrows():
        c = {"usgs": "usgs", "combined": "combined"}.get(r.catalog, "combined_nodup")
        for k in ["after_depth_filter", "M>=4.4", "M>=4.4_intermediate", "M>=4.4_shallow", "M>=4.4_uncertain",
                  "M>=4.4_deep", "all_intermediate", "all_shallow", "all_uncertain", "all_deep"]:
            if k in r:
                add(k, c, r.case, r[k])
    for c, base in CATS.items():
        for case in ("excluding_fixed", "including_fixed"):
            pp = read(os.path.join(base, case, "per_period_summary.csv"))
            if pp is not None:
                for _, r in pp.iterrows():
                    add(f"rate_per_yr_M4.4_{r.period}", c, case, round(r.rate_per_yr, 3))
                    add(f"b_MLE_{r.period}", c, case, round(r.b_MLE, 3))
            ks = read(os.path.join(base, case, "kishtwar", "kishtwar_summary.csv"))
            if ks is not None:
                for _, r in ks.iterrows():
                    add(f"kishtwar_{r.iloc[0]}", c, case, r.value)
            dc = read(os.path.join(base, case, "depth_clusters", "depth_clusters.csv"))
            if dc is not None:
                for g, n in dc.groupby("depth_group").size().items():
                    add(f"dbscan_clusters_{g}", c, case, n)
    bv = read(os.path.join(R, "robustness", "bvalues.csv"))
    for _, r in bv[bv.group.isin(["whole region", "intermediate"])].iterrows():
        add(f"b_mle_{r.group}_{r.magnitudes}", r.catalog, r.case, round(r.b_mle, 3))
        add(f"N_b_{r.group}_{r.magnitudes}", r.catalog, r.case, r.N)
    mc = read(os.path.join(R, "robustness", "mc_by_decade.csv"))
    for _, r in mc.iterrows():
        add(f"mc_authors_{r.decade}s", r.catalog, "", r.mc_authors)
    ms = read(os.path.join(R, "robustness", "mc_runs_summary.csv"))
    for _, r in ms.iterrows():
        add("mc_band_recovered_share_50pct", r.catalog, "excluding_fixed", r.band_recovered_share)
        add("mc_n_band_events", r.catalog, "excluding_fixed", r.get("n_band", r.get("n_in_outline")))
    bl = read(os.path.join(HERE, "results_baseline", "baseline_bvalues.csv"))
    for _, r in bl.iterrows():
        for k in bl.columns[1:4]:
            add(f"baseline_{r.group}_{k}", "authors", "", r[k])
        add(f"baseline_{r.group}_b_lsq", "authors", "", r.b_lsq)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"{len(rows)} metrics -> {out}")


if __name__ == "__main__":
    main(sys.argv[1])
