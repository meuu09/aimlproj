"""b-value robustness: magnitude types, error estimates and cutoff dependence (results/robustness/)."""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path

from seis_utils import (CLUSTER_POLY, HERE, b_bootstrap, b_lsq, b_mle, b_shi_bolt_sigma, load_catalog)

OUT = os.path.join(HERE, "results", "robustness")
START_YEAR, MC_CUT = 1990, 4.4
CUTOFFS = np.round(np.arange(4.2, 4.81, 0.1), 1)
os.makedirs(OUT, exist_ok=True)


def groups(df):
    """Event groups whose b-values are compared."""
    inside = Path(CLUSTER_POLY).contains_points(df[["longitude", "latitude"]].values)
    return {"whole region": df,
            "intermediate": df[df.depth_class == "intermediate"],
            "shallow": df[df.depth_class == "shallow"],
            "Kishtwar outline, all": df[inside],
            "Kishtwar outline, intermediate": df[inside & (df.depth_class == "intermediate")]}


def one(m, mc):
    """MLE b with Aki, Shi-Bolt and bootstrap errors, plus the LSQ b, above mc."""
    b, aki, n = b_mle(m, mc)
    med, lo, hi = b_bootstrap(m, mc)
    bl = b_lsq(m, mc)[0] if n >= 5 else np.nan
    return dict(N=n, b_mle=b, sigma_aki=aki, sigma_shi_bolt=b_shi_bolt_sigma(m, mc),
                boot_median=med, boot_ci95_low=lo, boot_ci95_high=hi, b_lsq=bl)


def main():
    """b-values per catalog, case, group and magnitude type, and b against the cutoff."""
    rows, curve = [], []
    for cat, src, nodup in (("usgs", "usgs", False), ("combined", "combined", False),
                            ("combined_nodup", "combined", True)):
        for inc in (False, True):
            df, _ = load_catalog(src, inc, START_YEAR, nodup, verbose=False)
            case = "including_fixed" if inc else "excluding_fixed"
            for gname, g in groups(df).items():
                for mags_used, sub in (("all types", g), ("mb only", g[g.magType_norm == "mb"])):
                    r = dict(catalog=cat, case=case, group=gname, magnitudes=mags_used, Mc=MC_CUT)
                    r.update(one(sub.mag.values, MC_CUT))
                    r["sample_size"] = ("N>=50" if r["N"] >= 50 else
                                        "N 20-49: rough" if r["N"] >= 20 else "N<20: too few")
                    rows.append(r)
                    if gname in ("whole region", "intermediate"):
                        for c in CUTOFFS:
                            b, aki, n = b_mle(sub.mag.values, c)
                            _, lo, hi = b_bootstrap(sub.mag.values, c, n_boot=1000)
                            curve.append(dict(catalog=cat, case=case, group=gname, magnitudes=mags_used,
                                              Mc=c, N=n, b_mle=b, boot_ci95_low=lo, boot_ci95_high=hi))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(OUT, "bvalues.csv"), index=False)
    cur = pd.DataFrame(curve)
    cur.to_csv(os.path.join(OUT, "bvalue_vs_cutoff.csv"), index=False)

    fig, axs = plt.subplots(1, 4, figsize=(17, 4.4), sharey=True)
    combos = [(c, k) for c in ("usgs", "combined") for k in ("excluding_fixed", "including_fixed")]
    style = {"whole region": dict(color="#1f3b73", marker="o"), "intermediate": dict(color="#c0392b", marker="s")}
    for ax, (c, k) in zip(axs, combos):
        for gname, st in style.items():
            d = cur[(cur.catalog == c) & (cur.case == k) & (cur.group == gname) & (cur.magnitudes == "mb only")]
            d = d[d.b_mle.notna()]
            off = -0.012 if gname == "whole region" else 0.012
            ax.errorbar(d.Mc + off, d.b_mle, yerr=[d.b_mle - d.boot_ci95_low, d.boot_ci95_high - d.b_mle],
                        fmt="-", lw=1.5, ms=5, capsize=2, elinewidth=1, label=gname, **st)
            for x, y, n in zip(d.Mc, d.b_mle, d.N):
                ax.annotate(str(n), (x + off, y), textcoords="offset points", xytext=(4, 4), fontsize=7,
                            color=st["color"])
        ax.axhline(1.0, color="0.6", lw=0.8, ls=":")
        ax.set_title(f"{'USGS' if c == 'usgs' else 'USGS + ISC'} | {k.replace('_', ' ')}", fontsize=10)
        ax.set_xlabel("Magnitude cutoff Mc (mb)"); ax.grid(alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
    axs[0].set_ylim(0.5, 4.0)
    axs[0].set_ylabel("b-value (MLE), bootstrap 95 % interval")
    axs[0].legend(fontsize=9, loc="upper left")
    fig.suptitle("b-value against magnitude cutoff, mb magnitudes only (numbers = events above the cutoff; "
                 "y-axis cut at 4, longer error bars clipped)",
                 fontsize=11)
    plt.tight_layout(); fig.savefig(os.path.join(OUT, "bvalue_vs_cutoff.png"), dpi=200); plt.close(fig)

    pd.set_option("display.width", 220)
    show = tab[tab.group.isin(["whole region", "intermediate", "Kishtwar outline, intermediate"])]
    print(show[["catalog", "case", "group", "magnitudes", "N", "b_mle", "sigma_aki", "sigma_shi_bolt",
                "boot_ci95_low", "boot_ci95_high", "b_lsq", "sample_size"]].round(2).to_string(index=False))
    print("Outputs in", OUT)


if __name__ == "__main__":
    main()
