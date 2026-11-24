"""Join metrics_before.csv and metrics_after.csv into BEFORE_AFTER_RESULTS.csv with an explanation per change."""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def num(x):
    """Float, or NaN for text."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def why(r):
    """Reason for a changed value, from the audit corrections that can cause it."""
    k, c = r.metric, r.catalog
    if pd.isna(r.value_before):
        return "new metric (added during the audit)"
    if pd.isna(r.value_after):
        return "metric renamed or removed (see CHANGELOG.md)"
    if not r.changed:
        return ""
    if k.startswith("baseline_"):
        return ("the old baseline CSV had unquoted commas, so its columns were shifted when read; it now reports the "
                "released-code Mc (Mc+0.4) separately from the smallest fitted magnitude (Task 1)")
    if k.startswith("mc_band") or k.startswith("mc_n_band"):
        return "joint depth+epicentre Monte Carlo, missing depth errors no longer set to zero (Tasks 18, 39, 40)"
    if k.startswith("b_MLE"):
        return "b-values from fewer than 20 events are no longer reported (Task 30)"
    if c.startswith("combined"):
        return ("combined rows now use one whole origin (Task 8), sentinel errors set missing, and the depth error is "
                "applied at both 40 and 150 km (Task 16)")
    if "rate" in k:
        return "catalog end date applied consistently (Task 13)"
    return "depth error applied at both 40 and 150 km bounds (Task 16)"


def main():
    """Write BEFORE_AFTER_RESULTS.csv."""
    b = pd.read_csv(os.path.join(HERE, "metrics_before.csv"))
    a = pd.read_csv(os.path.join(HERE, "metrics_after.csv"))
    for d in (a, b):
        d["case"] = d["case"].fillna("")
    m = b.merge(a, on=["metric", "catalog", "case"], how="outer", suffixes=("_before", "_after"))
    m["changed"] = [not (num(x) == num(y) or (np.isnan(num(x)) and np.isnan(num(y)) and str(x) == str(y)))
                    for x, y in zip(m.value_before, m.value_after)]
    m["explanation"] = m.apply(why, axis=1)
    m.sort_values(["changed", "metric", "catalog", "case"], ascending=[False, True, True, True]).to_csv(
        os.path.join(HERE, "BEFORE_AFTER_RESULTS.csv"), index=False)
    print(f"{len(m)} rows, {int(m.changed.sum())} changed")


if __name__ == "__main__":
    main()
