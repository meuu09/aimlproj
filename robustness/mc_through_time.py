"""Magnitude of completeness through time in sliding windows (results/robustness/)."""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from matplotlib.path import Path

from seis_utils import (CLUSTER_POLY, HERE, b_mle, b_reliability, load_catalog, max_curv, max_curv_eqcorrscan,
                        mc_gft, mc_mbs, mc_sliding)

OUT = os.path.join(HERE, "results", "robustness")
WINDOW, STEP, MC_CUT = 50, 10, 4.4
os.makedirs(OUT, exist_ok=True)


def mc_methods(us, co, nd):
    """Mc by five methods for each catalog, magnitude type and the Kishtwar outline, 1990 onward (Tasks 24-27)."""
    def inside(d):
        return d[Path(CLUSTER_POLY).contains_points(d[["longitude", "latitude"]].values)]
    groups = [("USGS", "all types", us), ("ISC-only rows", "all types", co[co.source == "ISC"]),
              ("USGS + ISC", "all types", co), ("USGS + ISC, duplicates removed", "all types", nd)]
    groups += [("USGS + ISC", t, co[co.magType_norm == t]) for t in ("mb", "ML", "MD", "Mw")]
    groups += [("USGS + ISC, duplicates removed", "mb", nd[nd.magType_norm == "mb"]),
               ("USGS, Kishtwar outline", "all types", inside(us)), ("USGS + ISC, Kishtwar outline", "mb",
                                                                     inside(co[co.magType_norm == "mb"]))]
    rows = []
    for name, mtype, d in groups:
        m = d[d.time.dt.year >= 1990].mag.values
        ok = len(m) >= 20
        auth = max_curv_eqcorrscan(m) + 0.4 if ok else np.nan
        frac = np.mean(np.abs(m * 10 - np.round(m * 10)) < 1e-6) if len(m) else np.nan
        b, s, n = b_mle(m, auth) if ok else (np.nan, np.nan, 0)
        rows.append(dict(catalog=name, magnitudes=mtype, N_1990on=len(m), share_on_0p1_grid=round(frac, 3),
                         mc_maxc=max_curv(m) if ok else np.nan, mc_maxc_plus0p4=max_curv(m) + 0.4 if ok else np.nan,
                         mc_authors=auth, mc_gft90=mc_gft(m, level=90), mc_gft95=mc_gft(m, level=95), mc_mbs=mc_mbs(m),
                         N_above_authors_mc=n, b_mle_at_authors_mc=b if n >= 20 else np.nan, b_reliability=b_reliability(n)))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, "mc_methods.csv"), index=False)
    print(t.round(2).to_string(index=False))


def main():
    """Mc in sliding windows and per decade for USGS and USGS + ISC."""
    us, _ = load_catalog("usgs", include_fixed=True, verbose=False)
    co, _ = load_catalog("combined", include_fixed=True, verbose=False)
    nd, _ = load_catalog("combined", include_fixed=True, drop_duplicates=True, verbose=False)
    sets = [("USGS (mostly mb)", us), ("USGS + ISC, mb only", co[co.magType_norm == "mb"]),
            ("USGS + ISC, all types (mixed scales)", co)]
    mc_methods(us, co, nd)
    rows, dec = [], []
    for name, df in sets:
        r = mc_sliding(df.time.values, df.mag.values, WINDOW, STEP)
        r.insert(0, "catalog", name)
        rows.append(r)
        for d, g in df.groupby(df.time.dt.year // 10 * 10):
            dec.append(dict(catalog=name, decade=int(d), N=len(g),
                            mc_maxc=max_curv(g.mag) if len(g) >= 20 else np.nan,
                            mc_authors=max_curv_eqcorrscan(g.mag.values) + 0.4 if len(g) >= 20 else np.nan,
                            N_M_ge_4p4=int((g.mag >= MC_CUT).sum())))
    tab = pd.concat(rows, ignore_index=True)
    tab.to_csv(os.path.join(OUT, "mc_through_time.csv"), index=False)
    dtab = pd.DataFrame(dec)
    dtab.to_csv(os.path.join(OUT, "mc_by_decade.csv"), index=False)

    fig, axs = plt.subplots(len(sets), 1, figsize=(12, 9), sharex=True)
    for ax, (name, df) in zip(axs, sets):
        d = tab[tab.catalog == name]
        ax.scatter(df.time, df.mag, s=4, c="0.75", lw=0, label="events")
        ax.fill_between(d.t_mid, d.mc_maxc_p16, d.mc_maxc_p84, color="#1f3b73", alpha=0.2, lw=0,
                        label="max. curvature, bootstrap 16-84 %")
        ax.plot(d.t_mid, d.mc_maxc, color="#1f3b73", lw=1.8, label="max. curvature (Wiemer & Wyss)")
        ax.plot(d.t_mid, d.mc_authors, color="#c0392b", lw=1.5, ls="--", label="authors' method (+0.4)")
        ax.axhline(MC_CUT, color="k", lw=0.8, ls=":")
        ax.set_title(f"{name}: windows of {WINDOW} events, step {STEP}", fontsize=10, loc="left")
        ax.set_ylabel("Magnitude"); ax.grid(alpha=0.2); ax.spines[["top", "right"]].set_visible(False)
    axs[0].legend(fontsize=8, loc="upper right", ncol=2)
    axs[-1].xaxis.set_major_locator(mdates.YearLocator(10)); axs[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    axs[-1].set_xlabel("Time (UTC); dotted line = Mc 4.4 used for rates and b-values")
    plt.tight_layout(); fig.savefig(os.path.join(OUT, "mc_through_time.png"), dpi=200); plt.close(fig)

    pd.set_option("display.width", 200)
    print(dtab.round(2).to_string(index=False))
    for name, _ in sets:
        d = tab[(tab.catalog == name) & (tab.t_mid.dt.year >= 1990)]
        print(f"{name}: windows from 1990, max-curvature Mc {d.mc_maxc.min():.1f}-{d.mc_maxc.max():.1f}, "
              f"authors' Mc {d.mc_authors.min():.1f}-{d.mc_authors.max():.1f}; "
              f"share of windows with authors' Mc > {MC_CUT}: {(d.mc_authors > MC_CUT + 1e-9).mean():.0%}")
    print("Outputs in", OUT)


if __name__ == "__main__":
    main()
