"""Catalog quality-control tables and QC_REPORT.md (results/qc/)."""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from seis_utils import (DUP_MAX_KM, DUP_MAX_S, HERE, SOLVED_DEPTH, close_pairs, load_catalog)

OUT = os.path.join(HERE, "results", "qc")
MC_CUT, START_YEAR = 4.4, 1990
REVIEW_S, REVIEW_KM = 60.0, 150.0
os.makedirs(OUT, exist_ok=True)


def decade(t):
    """Decade of each timestamp, e.g. 1990."""
    return (t.dt.year // 10 * 10).astype(int)


def main():
    """Write all QC tables, the match figure and QC_REPORT.md."""
    us, _ = load_catalog("usgs", include_fixed=True, verbose=False)
    co, _ = load_catalog("combined", include_fixed=True, verbose=False)
    lines = []

    mt_src = pd.crosstab(co.magType.fillna("none"), co.source, margins=True, margins_name="total")
    mt_src = mt_src.sort_values("total", ascending=False)
    mt_src.to_csv(os.path.join(OUT, "magtype_by_source.csv"))
    rows = []
    for name, df in [("usgs", us), ("combined", co)]:
        for thr in (None, MC_CUT):
            d = df if thr is None else df[df.mag >= thr]
            t = pd.crosstab(decade(d.time), d.magType.fillna("none"))
            t.insert(0, "subset", "all" if thr is None else f"M>={thr}")
            t.insert(0, "catalog", name)
            rows.append(t.reset_index().rename(columns={"time": "decade"}))
    pd.concat(rows, ignore_index=True).fillna(0).to_csv(os.path.join(OUT, "magtype_by_decade.csv"), index=False)
    c90 = co[co.time.dt.year >= START_YEAR]
    big = c90[c90.mag >= MC_CUT]
    lines += ["## Magnitude types", "",
              f"* Combined catalog, {START_YEAR}+: {len(c90)} events, magnitudes "
              f"{c90.mag.min():.1f} to {c90.mag.max():.1f}; {int((c90.mag < 3).sum())} are below M3.",
              "* Types (all years): " + ", ".join(f"{k} {v}" for k, v in co.magType.value_counts().items()) + ".",
              f"* At M >= {MC_CUT} ({START_YEAR}+): {len(big)} events, "
              f"{(big.magType == 'mb').mean():.0%} mb. Below about M4 the catalog is mostly ML/MD from "
              "India's NDI, so Mc and b-values computed below M4.4 mix scales and are not comparable.",
              "* No magnitude conversion is applied anywhere; magType and mag_source are kept.", ""]

    us["solved"] = us.depth_type.isin(SOLVED_DEPTH)
    dq = pd.crosstab(decade(us.time), us.depth_type)
    dq["N_M>=4.4"] = us[us.mag >= MC_CUT].groupby(decade(us[us.mag >= MC_CUT].time)).size()
    dq["share_solved_M>=4.4"] = us[us.mag >= MC_CUT].groupby(decade(us[us.mag >= MC_CUT].time)).solved.mean().round(2)
    dq.index.name = "decade"
    dq.to_csv(os.path.join(OUT, "depth_quality_by_decade.csv"))
    co["solved"] = co.depth_type.isin(SOLVED_DEPTH)
    ds = pd.crosstab(decade(co[co.solved].time), co[co.solved].depth_source)
    ds.index.name = "decade"
    ds.to_csv(os.path.join(OUT, "depth_source_by_decade.csv"))
    lines += ["## Depth quality through time (USGS)", "",
              "| decade | events M>=4.4 | share with solved depth |", "|---|---|---|"]
    lines += [f"| {d} | {int(r['N_M>=4.4'])} | {r['share_solved_M>=4.4']:.0%} |" for d, r in dq.iterrows()]
    lines += ["", "The share of events whose depth was actually solved changes a lot between decades. "
              "In the excluding-fixed case, an event's chance of being kept therefore depends on its decade, "
              "so event rates per period in that case partly measure USGS practice, not the Earth. "
              "Do not interpret rate changes between periods as physical.", "",
              "Combined catalog, solved depths by agency (all years): "
              + ", ".join(f"{k} {v}" for k, v in co[co.solved].depth_source.value_counts().items())
              + ". For an ISC depth, the epicentre and depth may come from different agencies "
              "(epicentre: USGS or ISC prime; depth: first trusted agency with a free depth and an error).", ""]

    fun = []
    for cat_name, dd in [("usgs", False), ("combined", False), ("combined", True)]:
        for inc in (False, True):
            df, info = load_catalog(cat_name, inc, START_YEAR, dd, verbose=False)
            big = df[df.mag >= MC_CUT]
            row = dict(catalog=cat_name + (" (duplicates removed)" if dd else ""),
                       case="including_fixed" if inc else "excluding_fixed", **{k: v for k, v in info.items()
                                                                                if k != "catalog"})
            row["M>=4.4"] = len(big)
            for k in ("intermediate", "shallow", "uncertain", "deep"):
                row[f"M>=4.4_{k}"] = int((big.depth_class == k).sum())
                row[f"all_{k}"] = int((df.depth_class == k).sum())
            fun.append(row)
    fun = pd.DataFrame(fun)
    fun.to_csv(os.path.join(OUT, "filter_funnel.csv"), index=False)
    lines += ["## Filtering summary (events left after each step)", "",
              "| catalog | case | in box | " + f"{START_YEAR}+ | after dedup | after depth filter | M>=4.4 | "
              "of which intermediate / shallow / uncertain / >150 km |", "|---|---|---|---|---|---|---|---|"]
    for _, r in fun.iterrows():
        lines.append(f"| {r.catalog} | {r.case} | {r.in_box_with_mag} | {r.after_start_year} | {r.after_duplicates} | "
                     f"{r.after_depth_filter} | {r['M>=4.4']} | {r['M>=4.4_intermediate']} / {r['M>=4.4_shallow']} / "
                     f"{r['M>=4.4_uncertain']} / {r['M>=4.4_deep']} |")
    lines.append("")

    m = co[co.source == "USGS+ISC"].copy()
    m["dmag"] = m.mag - m.isc_mag
    m["same_magType"] = m.magType.str.lower() == m.isc_magType.str.lower()
    cols = ["id", "isc_id", "time", "match_dt_s", "match_dist_km", "mag", "magType", "isc_mag", "isc_magType",
            "dmag", "same_magType"]
    m[cols].to_csv(os.path.join(OUT, "match_diagnostics.csv"), index=False)
    q = m[["match_dt_s", "match_dist_km", "dmag"]].describe(percentiles=[0.5, 0.95]).round(2)
    pre, post = m[m.time.dt.year < 2000], m[m.time.dt.year >= 2000]
    lines += ["## USGS-ISC matching (window 20 s, 100 km, one-to-one)", "",
              f"* {len(m)} matches. Time difference median {q.loc['50%', 'match_dt_s']} s, 95th percentile "
              f"{q.loc['95%', 'match_dt_s']} s, max {q.loc['max', 'match_dt_s']} s.",
              f"* Distance median {q.loc['50%', 'match_dist_km']} km, 95th percentile {q.loc['95%', 'match_dist_km']} km, "
              f"max {q.loc['max', 'match_dist_km']} km; {int((m.match_dist_km > 50).sum())} matches are more than 50 km "
              f"apart (median distance before 2000: {pre.match_dist_km.median():.1f} km, from 2000: "
              f"{post.match_dist_km.median():.1f} km).",
              f"* Magnitude difference (USGS minus ISC-chosen) median {q.loc['50%', 'dmag']}; "
              f"for matches with the same magnitude type: median "
              f"{m[m.same_magType].dmag.median():.2f}, |difference| > 0.5 in "
              f"{int((m[m.same_magType].dmag.abs() > 0.5).sum())} of {int(m.same_magType.sum())}.",
              "* The magnitudes agree closely, which supports the matches being the same earthquakes. "
              "Matches more than 50 km apart by decade: "
              + ", ".join(f"{k}s {v}" for k, v in decade(m[m.match_dist_km > 50].time).value_counts().sort_index().items())
              + " (listed in match_diagnostics.csv; check them by hand before relying on those events).", ""]
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, col, lab, bins in [(axs[0], "match_dt_s", "Origin-time difference (s)", np.arange(0, 21, 1)),
                               (axs[1], "match_dist_km", "Epicentre distance (km)", np.arange(0, 105, 5)),
                               (axs[2], "dmag", "Magnitude difference, USGS − ISC", np.arange(-1.05, 1.1, 0.1))]:
        ax.hist(m[col].dropna(), bins=bins, color="#4a6fa5", edgecolor="white", linewidth=1)
        ax.set_xlabel(lab); ax.grid(axis="y", alpha=0.25); ax.spines[["top", "right"]].set_visible(False)
    axs[0].set_ylabel("Matched events")
    fig.suptitle(f"USGS–ISC matches (N = {len(m)}); match window 20 s and 100 km", fontsize=11)
    plt.tight_layout(); fig.savefig(os.path.join(OUT, "match_diagnostics.png"), dpi=200); plt.close(fig)

    pr = close_pairs(co, REVIEW_S, REVIEW_KM)
    a, b = co.iloc[pr.i.to_numpy(int)].reset_index(drop=True), co.iloc[pr.j.to_numpy(int)].reset_index(drop=True)
    cp = pd.DataFrame({"time_1": a.time, "id_1": a.id, "source_1": a.source, "mag_1": a.mag, "magType_1": a.magType,
                       "mag_source_1": a.mag_source, "id_2": b.id, "source_2": b.source, "mag_2": b.mag,
                       "magType_2": b.magType, "mag_source_2": b.mag_source,
                       "dt_s": pr.dt_s.round(2), "dist_km": pr.dist_km.round(1)})
    cp["flagged"] = (a.possible_duplicate_of.notna() | b.possible_duplicate_of.notna()).values
    cp["why"] = np.where(cp.flagged, f"within {DUP_MAX_S:.0f} s and {DUP_MAX_KM:.0f} km, at least one ISC-only",
                         np.where((cp.source_1 != "ISC") & (cp.source_2 != "ISC"),
                                  "both from USGS (listed separately by USGS): not flagged",
                                  f"farther than {DUP_MAX_KM:.0f} km: not flagged"))
    cp.to_csv(os.path.join(OUT, "close_pairs.csv"), index=False)
    nflag = int(co.possible_duplicate_of.notna().sum())
    fl = cp[cp.flagged]
    lines += ["## Possible duplicate events", "",
              f"* {len(cp)} pairs of rows lie within {REVIEW_S:.0f} s and {REVIEW_KM:.0f} km of each other "
              f"(close_pairs.csv). {nflag} rows are flagged as possible duplicates "
              f"(rule: within {DUP_MAX_S:.0f} s and {DUP_MAX_KM:.0f} km, at least one row ISC-only; "
              "column possible_duplicate_of).",
              f"* Flagged pairs: origin-time differences {fl.dt_s.min():.1f}-{fl.dt_s.max():.1f} s, distances "
              f"{fl.dist_km.min():.1f}-{fl.dist_km.max():.1f} km; magnitude sources of the flagged rows: "
              + ", ".join(f"{k} {v}" for k, v in co[co.possible_duplicate_of.notna()].mag_source.value_counts().items())
              + ".",
              "* Flagged rows stay in the catalog. The analysis scripts leave them out only with "
              "--drop-duplicates (results in results/combined_usgs_isc_nodup/).", ""]

    head = ["# Catalog quality report", "",
            "Generated by catalog_qc.py from data/ (no data changed). Tables are in this folder.", ""]
    open(os.path.join(OUT, "QC_REPORT.md"), "w", encoding="utf-8").write("\n".join(head + lines) + "\n")
    print("\n".join(lines))
    print("Outputs in", OUT)


if __name__ == "__main__":
    main()
