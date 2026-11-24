"""
Magnitude of completeness through time (results/robustness/).

Rates can only be compared between periods above a magnitude that was complete in
ALL of them. This script estimates Mc in sliding windows of consecutive events
(50 events, moved 10 at a time), with two methods:
  * maximum curvature (Wiemer & Wyss 2000) with a bootstrap 16-84 % range
  * the authors' eqcorrscan 0.4.2 method + 0.4 (as in Michailos et al. 2021)

Completeness is about detection, not depth quality, so all events are used
(fixed depths included), all years. Catalogs:
  * USGS (M >= 3.1 download, so Mc cannot come out below ~3.1)
  * USGS + ISC, mb only (one magnitude scale)
  * USGS + ISC, all types: shown for completeness only, it mixes mb, ML, MD and Mw

Outputs: mc_through_time.csv, mc_through_time.png, mc_by_decade.csv
Run:  .venv/bin/python mc_through_time.py      Windows: .venv\\Scripts\\python mc_through_time.py
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from seis_utils import HERE, load_catalog, max_curv, max_curv_eqcorrscan, mc_sliding

OUT = os.path.join(HERE, "results", "robustness")
WINDOW, STEP, MC_CUT = 50, 10, 4.4
os.makedirs(OUT, exist_ok=True)


def main():
    us, _ = load_catalog("usgs", include_fixed=True, verbose=False)
    co, _ = load_catalog("combined", include_fixed=True, verbose=False)
    sets = [("USGS (mostly mb)", us), ("USGS + ISC, mb only", co[co.magType == "mb"]),
            ("USGS + ISC, all types (mixed scales)", co)]
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
