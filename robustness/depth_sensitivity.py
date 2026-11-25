"""How counts, clusters and temporal statistics change with the intermediate-depth bounds (results/robustness/)."""
import os

import numpy as np
import pandas as pd
from matplotlib.path import Path
from sklearn.cluster import DBSCAN

from seis_utils import CATALOG_END, CLUSTER_POLY, EARTH_R, HERE, cv_poisson_test, depth_class, load_catalog

OUT = os.path.join(HERE, "results", "robustness")
LOWER, UPPER = [40, 50, 60, 70], [100, 150]
os.makedirs(OUT, exist_ok=True)


def main():
    """Reclassify with every (lower, upper) pair and summarise."""
    rows = []
    for name, src, nodup in (("usgs", "usgs", False), ("combined", "combined", False),
                             ("combined_nodup", "combined", True)):
        for inc in (False, True):
            df, _ = load_catalog(src, inc, None, nodup, verbose=False)
            inside = Path(CLUSTER_POLY).contains_points(df[["longitude", "latitude"]].values)
            for lo in LOWER:
                for up in UPPER:
                    cls = depth_class(df.depth, df.depthError, lo, up)
                    inter = cls == "intermediate"
                    d = df[inter]
                    lab = (DBSCAN(eps=30 / EARTH_R, min_samples=4, metric="haversine")
                           .fit_predict(np.radians(d[["latitude", "longitude"]].values)) if len(d) >= 4
                           else np.full(len(d), -1))
                    recent = d[d.time.dt.year >= 1990]
                    c, c_lo, c_hi, p = cv_poisson_test(recent.time, "1990-01-01", CATALOG_END)
                    rows.append(dict(catalog=name, case="including_fixed" if inc else "excluding_fixed",
                                     lower_km=lo, upper_km=up, n_intermediate=int(inter.sum()),
                                     n_uncertain=int((cls == "uncertain").sum()),
                                     n_intermediate_M4p4_1990on=int(((d.mag >= 4.4) & (d.time.dt.year >= 1990)).sum()),
                                     n_in_kishtwar_outline=int((inter & inside).sum()),
                                     n_dbscan_clusters=len(set(lab[lab >= 0])),
                                     largest_cluster=int(np.bincount(lab[lab >= 0]).max()) if (lab >= 0).any() else 0,
                                     depth_median=d.depth.median(), cv_1990on=c, cv_poisson_p=p))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, "depth_threshold_sensitivity.csv"), index=False)
    pd.set_option("display.width", 200)
    print(t[t.case == "excluding_fixed"].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
