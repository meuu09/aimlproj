"""Compare the reproduced Him_ides results with the values published in Michailos et al. (2021)."""
import os

import numpy as np
import pandas as pd

from seis_utils import HERE, b_lsq, b_mle, max_curv_eqcorrscan

OUT = os.path.join(HERE, "results_baseline")
SRC = "Michailos et al. 2021, Sect. 3.2.2"


def status(diff, tol):
    """Exact, approximate or unresolved, from an absolute difference and a tolerance."""
    if np.isnan(diff):
        return "not comparable"
    return "exact" if abs(diff) < 1e-9 else ("approximate" if abs(diff) <= tol else "unresolved discrepancy")


def main():
    """Write baseline_comparison.csv."""
    b = pd.read_csv(os.path.join(OUT, "baseline_bvalues.csv")).set_index("group")
    mags = pd.read_csv(os.path.join(OUT, "baseline_magnitudes.csv"))
    full = mags[mags.group == "Full catalog"].mag.to_numpy()
    c1 = mags[mags.group == "Cluster 1"].mag.to_numpy()
    rows = [
        ("Full catalog: number of events", 414, len(full), 0, "Abstract and Sect. 3.2.2",
         "Same bundled catalog."),
        ("Full catalog: smallest ML", 0.79, full.min(), 0.005, SRC,
         "One bundled event has ML 0.14; the next smallest is 0.79, as published. The released file holds one event "
         "(or magnitude) not in the published range."),
        ("Full catalog: largest ML", 4.51, full.max(), 0.005, SRC, ""),
        ("Full catalog: Mc (max curvature + 0.4)", 2.4, b.loc["Full catalog", "Mc_released_code"], 0.0, SRC,
         "Released calc_max_curv (eqcorrscan 0.4.2, 0.5 bins) + 0.4."),
        ("Full catalog: b (fit to cumulative distribution)", 0.91, b.loc["Full catalog", "b_lsq"], 0.004, SRC,
         "Released least-squares fit; 0.8 % above the published value, outside its stated +/-0.004. That error is a "
         "formal fit error, far smaller than the sampling error "
         f"(Aki b/sqrt(N) = {b_mle(full, 2.4, dm=0.01)[1]:.3f})."),
        ("Cluster 1: number of events", 248, b.loc["Cluster 1", "N"], 0, "Sect. 3.1",
         "Paper selects 'box in Figure 2'; released code uses two polygons (NE of Everest 1 + SE of one), which "
         "give 254 events whichever origin (first, last, preferred) is used. The exact box is not in the code."),
        ("Cluster 1: Mc", 2.85, b.loc["Cluster 1", "Mc_released_code"], 0.0, SRC,
         "Released code bins at 0.5, so Mc + 0.4 can only end in .4 or .9; 2.85 cannot come from it. Bins of "
         f"0.2 / 0.1 / 0.05 give {max_curv_eqcorrscan(c1, 0.2) + 0.4:.2f} / {max_curv_eqcorrscan(c1, 0.1) + 0.4:.2f} / "
         f"{max_curv_eqcorrscan(c1, 0.05) + 0.4:.2f}. Not tuned to match."),
        ("Cluster 1: b at released Mc (2.4)", 0.96, b.loc["Cluster 1", "b_lsq"], 0.02, SRC,
         "Released LSQ fit at the released Mc."),
        ("Cluster 1: b, released LSQ fit at published Mc 2.85", 0.96, b_lsq(c1, 2.85)[0], 0.02, SRC,
         "Same fit, threshold set to the published Mc: shows the b difference comes from the threshold."),
        ("Cluster 1: b, maximum likelihood at Mc 2.85 (our check)", 0.96, b_mle(c1, 2.85, dm=0.01)[0],
         b_mle(c1, 2.85, dm=0.01)[1], SRC, f"Aki-Utsu MLE, N = {b_mle(c1, 2.85, dm=0.01)[2]}; tolerance = 1 sigma."),
        ("Cluster 2: number of events", 61, b.loc["Cluster 2", "N"], 0, SRC, "Udayapur polygon."),
        ("Cluster 2: b", np.nan, b.loc["Cluster 2", "b_lsq"], np.nan, SRC,
         "Paper: too few events for a robust b; the released code still plots one."),
        ("Migration sequences (20 d, 20 km, 3-D)", 4, np.nan, np.nan, "Sect. 3.2.3",
         "The search code is not in the released repository (supplementary material only); not reproduced."),
    ]
    t = pd.DataFrame(rows, columns=["parameter", "published", "reproduced", "tolerance", "source", "explanation"])
    t["abs_difference"] = (t.reproduced - t.published).abs().round(3)
    t["pct_difference"] = (100 * (t.reproduced - t.published) / t.published).round(1)
    t["status"] = [status(r.reproduced - r.published, r.tolerance) for r in t.itertuples()]
    t["reproduced"] = t.reproduced.round(3)
    t = t[["parameter", "published", "reproduced", "abs_difference", "pct_difference", "tolerance", "status",
           "source", "explanation"]]
    t.to_csv(os.path.join(OUT, "baseline_comparison.csv"), index=False)
    print(t[["parameter", "published", "reproduced", "status"]].to_string(index=False))


if __name__ == "__main__":
    main()
