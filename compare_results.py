"""
Compare the summary tables in results/ with the same files at an earlier git commit,
so every changed number is visible. Writes results/qc/old_vs_new.csv.

  .venv/bin/python compare_results.py            (compares with HEAD, i.e. the last commit)
  .venv/bin/python compare_results.py <commit>   (any other commit or branch)

Files compared (wherever they exist in both versions):
  per_period_summary.csv, kishtwar_summary.csv, kishtwar_bvalues.csv, depth_clusters.csv,
  and the event count / depth-class counts of catalog_M4.4_1990plus.csv.
Columns renamed on purpose are mapped (N_M3 -> N_all_used, diffusivity column).
"""
import io
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REF = sys.argv[1] if len(sys.argv) > 1 else "HEAD"
OUT = os.path.join(HERE, "results", "qc", "old_vs_new.csv")
KEYS = {"per_period_summary.csv": ["period"], "kishtwar_bvalues.csv": ["group"],
        "depth_clusters.csv": ["depth_group", "cluster"], "kishtwar_summary.csv": None}
RENAMED = {"N_M3": "N_all_used", "diffusivity_envelope_m2_per_s": "diffusivity_min_envelope_m2_per_s"}


def git_files(ref):
    out = subprocess.run(["git", "ls-tree", "-r", "--name-only", ref, "results"], cwd=HERE,
                         capture_output=True, text=True, check=True).stdout
    return out.split()


def git_read(ref, path):
    txt = subprocess.run(["git", "show", f"{ref}:{path}"], cwd=HERE, capture_output=True, text=True,
                         check=True, encoding="utf-8").stdout
    return pd.read_csv(io.StringIO(txt))


def long_form(df, name):
    if name == "kishtwar_summary.csv":
        df = df.rename(columns={df.columns[0]: "key"}).set_index("key")[["value"]].T
        df.index = ["summary"]
    else:
        df = df.rename(columns=RENAMED)
        keys = KEYS[name]
        if not set(keys) <= set(df.columns):
            return pd.DataFrame(columns=["row", "column", "value"])
        df = df.set_index(keys)
        df.index = [" / ".join(map(str, k)) if isinstance(k, tuple) else str(k) for k in df.index]
    df = df.rename(columns=RENAMED)
    s = df.stack(future_stack=True)
    return pd.DataFrame({"row": s.index.get_level_values(0), "column": s.index.get_level_values(1),
                         "value": s.values})


def catalog_counts(df):
    d = {"n_events": len(df)}
    cls = df.depth_class.astype(str).str.split(" ").str[0].value_counts()
    d.update({f"n_{k}": v for k, v in cls.items()})
    return pd.DataFrame({"row": "catalog", "column": list(d), "value": list(d.values())})


def same(a, b):
    try:
        fa, fb = float(a), float(b)
        return (np.isnan(fa) and np.isnan(fb)) or np.isclose(fa, fb, rtol=1e-6, atol=1e-9)
    except (TypeError, ValueError):
        return str(a) == str(b)


def main():
    rows = []
    for path in git_files(REF):
        name = os.path.basename(path)
        if name not in KEYS and not name.startswith("catalog_M"):
            continue
        new_path = os.path.join(HERE, path)
        old = git_read(REF, path)
        new = pd.read_csv(new_path) if os.path.exists(new_path) else None
        if new is None:
            rows.append(dict(file=path, row="", column="", old="(file)", new="(missing now)", changed=True))
            continue
        o, n = ((catalog_counts(old), catalog_counts(new)) if name.startswith("catalog_M")
                else (long_form(old, name), long_form(new, name)))
        m = o.merge(n, on=["row", "column"], how="outer", suffixes=("_old", "_new"))
        for r in m.itertuples(index=False):
            rows.append(dict(file=path, row=r.row, column=r.column, old=r.value_old, new=r.value_new,
                             changed=not same(r.value_old, r.value_new)))
    tab = pd.DataFrame(rows)
    ref_hash = subprocess.run(["git", "rev-parse", "--short", REF], cwd=HERE, capture_output=True, text=True,
                              check=True).stdout.strip()
    tab.insert(0, "compared_with", ref_hash)
    tab.to_csv(OUT, index=False)
    ch = tab[tab.changed]
    print(f"Compared {tab.file.nunique()} files with {REF}: {len(tab)} values, {len(ch)} changed -> {OUT}")
    pd.set_option("display.width", 220); pd.set_option("display.max_colwidth", 60)
    print(ch.to_string(index=False) if len(ch) else "no changes")


if __name__ == "__main__":
    main()
