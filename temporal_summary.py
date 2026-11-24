"""Temporal summary per depth population and completeness scenario, with rate CIs and Poisson tests (results/temporal/)."""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from seis_utils import (CATALOG_END, HERE, SOLVED_DEPTH, b_mle, b_reliability, cv_poisson_test, dispersion_index,
                        load_catalog, max_curv_eqcorrscan, poisson_ci)

OUT = os.path.join(HERE, "results", "temporal")
CATALOGS = [("usgs", "usgs", False), ("combined", "combined", False), ("combined_nodup", "combined", True)]
SCENARIOS = [("M4.4 from 1990 (exploratory)", 4.4, 1990, None),
             ("M4.9 from 1990 (conservative)", 4.9, 1990, None),
             ("M4.4 from 2000", 4.4, 2000, None),
             ("mb 3.9 from 1990 (mb only)", 3.9, 1990, "mb")]
POPULATIONS = ["all", "intermediate", "shallow", "uncertain", "deep"]
DECADES = [(1990, 1999), (2000, 2009), (2010, 2019), (2020, 2026)]
os.makedirs(OUT, exist_ok=True)


def decade_mc(cat_name, df, mag_type):
    """Authors' Mc (+0.4) per decade, from all events of that catalog (fixed depths included)."""
    d = df if mag_type is None else df[df.magType_norm == mag_type]
    return {int(k): max_curv_eqcorrscan(g.mag.values) + 0.4 for k, g in d.groupby(d.time.dt.year // 10 * 10)
            if len(g) >= 20}


def summarise(df, start, end, mc):
    """Count, duration, rate with exact 95 % CI, CV test, dispersion and b (when N allows) for one population."""
    years = (end - start).days / 365.25
    n = len(df)
    lo, hi = poisson_ci(n, years)
    c, c_lo, c_hi, p = cv_poisson_test(df.time, start, end)
    d, d_p = dispersion_index(df.time, start, end - pd.Timedelta(seconds=1))
    b, s, nb = b_mle(df.mag.values, mc)
    rel = b_reliability(nb)
    if nb < 20:
        b = s = np.nan
    return dict(N=n, years=round(years, 2), rate_per_yr=n / years, rate_ci95_low=lo, rate_ci95_high=hi,
                cv=c, cv_poisson_2p5=c_lo, cv_poisson_97p5=c_hi, cv_p_two_sided=p,
                annual_dispersion=d, dispersion_p=d_p, b_mle=b, b_sigma_aki=s, b_N=nb, b_reliability=rel)


def main():
    """Write temporal_summary.csv, intermediate_summary.csv and the rate figure."""
    rows = []
    for name, cat, nodup in CATALOGS:
        full, _ = load_catalog(cat, True, None, nodup, verbose=False)
        for case in ("excluding_fixed", "including_fixed"):
            df = full if case == "including_fixed" else full[full.depth_type.isin(SOLVED_DEPTH)]
            for scen, thr, y0, mtype in SCENARIOS:
                if mtype and name == "usgs":
                    continue
                mcs = decade_mc(name, full, mtype)
                sel = df[(df.mag >= thr - 1e-9) & (df.time.dt.year >= y0)]
                ref = full[(full.mag >= thr - 1e-9) & (full.time.dt.year >= y0)]
                if mtype:
                    sel, ref = sel[sel.magType_norm == mtype], ref[ref.magType_norm == mtype]
                for a, b in [(y0, 2026)] + [p for p in DECADES if p[0] >= y0]:
                    start = pd.Timestamp(f"{a}-01-01")
                    end = min(pd.Timestamp(f"{b + 1}-01-01"), CATALOG_END)
                    seg = sel[(sel.time >= start) & (sel.time < end)]
                    rseg = ref[(ref.time >= start) & (ref.time < end)]
                    worst = max([v for k, v in mcs.items() if a <= k + 9 and k <= b] or [np.nan])
                    for pop in POPULATIONS:
                        part = seg if pop == "all" else seg[seg.depth_class == pop]
                        r = dict(catalog=name, case=case, scenario=scen, threshold=thr, magnitudes=mtype or "all types",
                                 period=f"{a}-{b}", population=pop, max_decade_mc_authors=round(worst, 2),
                                 complete_by_authors_mc=bool(worst <= thr + 1e-9),
                                 share_solved_depth=(rseg.depth_type.isin(SOLVED_DEPTH).mean() if len(rseg) else np.nan))
                        r.update(summarise(part, start, end, thr))
                        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, "temporal_summary.csv"), index=False)
    inter = t[t.population == "intermediate"]
    inter.to_csv(os.path.join(OUT, "intermediate_summary.csv"), index=False)
    plot(t)
    pd.set_option("display.width", 250)
    cols = ["catalog", "case", "scenario", "period", "N", "rate_per_yr", "rate_ci95_low", "rate_ci95_high", "cv",
            "cv_p_two_sided", "complete_by_authors_mc", "share_solved_depth"]
    show = inter[(inter.case == "excluding_fixed") & (inter.catalog != "combined_nodup")]
    print(show[cols].round(2).to_string(index=False))
    print("Outputs in", OUT)


def plot(t):
    """Intermediate vs all-event rates per decade with 95 % CIs, one panel per catalog and scenario."""
    sub = t[(t.case == "excluding_fixed") & (t.catalog != "combined_nodup") & (t.period.str[:4] != t.period.str[5:])]
    sub = sub[sub.period.isin([f"{a}-{b}" for a, b in DECADES])]
    panels = sub[["catalog", "scenario"]].drop_duplicates().values.tolist()
    fig, axs = plt.subplots(1, len(panels), figsize=(4.2 * len(panels), 4.2), sharey=False)
    for ax, (cat, scen) in zip(np.atleast_1d(axs), panels):
        for pop, col, off in [("all", "#555555", -0.1), ("intermediate", "#c62828", 0.1)]:
            d = sub[(sub.catalog == cat) & (sub.scenario == scen) & (sub.population == pop)]
            x = np.arange(len(d)) + off
            ax.errorbar(x, d.rate_per_yr, yerr=[d.rate_per_yr - d.rate_ci95_low, d.rate_ci95_high - d.rate_per_yr],
                        fmt="o", color=col, capsize=3, label=pop)
            for xi, (_, r) in zip(x, d.iterrows()):
                if not r.complete_by_authors_mc:
                    ax.annotate("inc.", (xi, r.rate_ci95_high), fontsize=7, ha="center", va="bottom", color="0.4")
        ax.set_xticks(np.arange(len(d)), d.period, rotation=25, fontsize=8)
        ax.set_title(f"{'USGS' if cat == 'usgs' else 'USGS + ISC'}\n{scen}", fontsize=9)
        ax.set_yscale("log"); ax.grid(alpha=0.2)
    np.atleast_1d(axs)[0].set_ylabel("Events per year (exact Poisson 95 % CI)")
    np.atleast_1d(axs)[0].legend(fontsize=8)
    fig.suptitle("Solved depths only. 'inc.' = threshold below the authors' Mc for some decade in that period "
                 "(rate not comparable)", fontsize=9)
    plt.tight_layout(); fig.savefig(os.path.join(OUT, "rates_by_decade.png"), dpi=200); plt.close(fig)


if __name__ == "__main__":
    main()
