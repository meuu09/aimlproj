# Temporal evolution of intermediate-depth seismicity beneath the NW (Indian) Himalaya

IIT Roorkee group assignment (topic 4). We adapt the code of
**Michailos et al. (2021)**, *Frontiers in Earth Science* ([kemichai/Him_ides](https://github.com/kemichai/Him_ides),
vendored unchanged in `Him_ides/` at commit 33acd38; the repository has no licence file, so it is
credited here and in the report), reproduce their results, and apply the method to the
Himachal–Kashmir Himalaya (31–35°N, 74–79.5°E) using open data from **USGS ComCat** and the
**ISC Bulletin**.

## Setup

Python 3.12 is required. The original code imports eqcorrscan, which does not install on 3.12,
so `reproduce_baseline.py` supplies its one needed function (copied verbatim) instead.

**Windows (PowerShell)**

```powershell
py -3.12 -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m pytest -q tests
.venv\Scripts\python run_all.py
```

**macOS / Linux**

```bash
python3.12 -m venv .venv          # or: uv venv --python 3.12 .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pytest -q tests
.venv/bin/python run_all.py
```

`run_all.py` rebuilds every result from the files already in `data/` (under a minute, no downloads;
add `--baseline` to also rerun the original Him_ides code). Results are deterministic: two runs give
byte-identical CSV files.
The downloads themselves (`fetch_depth_type.py`, `fetch_isc_bulletin.py`) only need re-running to refresh the data.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 1. USGS events (M≥3, downloaded once from the FDSN API) | – | `data/usgs_himachal_kashmir_M3plus_1970-2026.csv` |
| 2. USGS depth flag per event (solved vs fixed by an analyst) | `fetch_depth_type.py` | `data/usgs_depth_type.csv` |
| 3. ISC Bulletin, all agencies' origins (ISF text) | `fetch_isc_bulletin.py` | `data/isc_raw/`, `data/isc_bulletin_events.csv`, `_origins.csv` |
| 4. Merge USGS + ISC, flag possible duplicates | `build_combined_catalog.py` | `data/combined_usgs_isc_catalog.csv` |
| 5. Catalog quality reports | `catalog_qc.py` | `results/qc/` (start with `QC_REPORT.md`) |
| 6. Regional temporal analysis (Figs 1–6) | `him_ides_himachal.py` | `results/<catalog>/<case>/` |
| 7. Kishtwar–Chamba cluster (paper's Figs 3–9 analogues) | `kishtwar_cluster.py` | `.../<case>/kishtwar/` |
| 8. DBSCAN clusters per depth class | `depth_clusters.py` | `.../<case>/depth_clusters/` |
| 9. b-value robustness (mb only, Shi & Bolt, bootstrap, cutoff) | `bvalue_robustness.py` | `results/robustness/` |
| 10. Mc through time (sliding windows) | `mc_through_time.py` | `results/robustness/` |
| 11. Cluster robustness (DBSCAN grid, Monte Carlo over errors) | `cluster_robustness.py` | `results/robustness/` |
| Baseline: original code on the authors' catalog (Figs 5–10) | `reproduce_baseline.py --fig10` | `results_baseline/` |
| Old vs new numbers after a code change | `compare_results.py [commit]` | `results/qc/old_vs_new.csv` |

Shared functions (loading, depth classes, Mc, b-value, matching, duplicate flags) are in
`seis_utils.py`; tests are in `tests/`.

### Options of the analysis scripts (steps 6–8)

| Option | Meaning | Results folder |
|---|---|---|
| *(none)* | USGS, **excluding** fixed/unflagged depths | `results/excluding_fixed/` |
| `--include-fixed` | USGS, all depths as reported | `results/including_fixed/` |
| `--catalog combined` | USGS + ISC Bulletin | `results/combined_usgs_isc/<case>/` |
| `--catalog combined --drop-duplicates` | same, without rows flagged as possible duplicates | `results/combined_usgs_isc_nodup/<case>/` |

## Definitions used everywhere

* **Depth classes:** intermediate = depth − depthError ≥ 40 km and depth ≤ 150 km; shallow = depth + depthError < 40 km;
  deep = depth > 150 km; uncertain = everything else (error bar crosses 40 km, or no error).
* **Solved depth:** USGS depth-type "from location" or "constrained by depth phases"; for ISC rows, a free
  (not fixed) depth with an error from ISC-EHB, ISC, NEIC or NDI.
* **Magnitudes are not converted.** USGS is mostly mb; ISC-only events are mostly ML/MD (India's NDI).
  `magType` and `mag_source` are kept. Rates and b-values use M ≥ 4.4, where 89–91% of events are mb.
* **Mc:** maximum curvature, two versions: Wiemer & Wyss (2000) mode, and the authors' eqcorrscan 0.4.2
  method + 0.4. b-value: Aki (1965)/Utsu maximum likelihood (preferred) and the repo's least-squares fit.

## What the results support (details and numbers in the files named)

* **Supported:** a band of intermediate-depth earthquakes (median depth ~60 km) beneath the
  Kishtwar–Chamba area. The hand-drawn outline and DBSCAN find it in both catalogs. At DBSCAN
  radii of 30–45 km, 64–100 % of its events fall in one cluster for every minimum size tested
  (USGS alone at 30 km splits it into a Kishtwar and a Chamba part), and at least half of them
  stay in one cluster in 99–100 % of Monte Carlo runs with epicentres shifted within their errors
  (`results/robustness/dbscan_grid.csv`, `mc_runs_summary.csv`). At a 15 km radius it falls
  apart, so its shape is only resolved at ~20–30 km, about 2–3 times the location errors.
* **Supported:** the original paper's numbers are reproduced (`results_baseline/`).
* **Not supported:** changes in earthquake rate over time. USGS alone is not complete at M4.4
  before ~2000 (`results/robustness/mc_by_decade.csv`), and the share of solved depths changes by
  decade (`results/qc/QC_REPORT.md`).
* **Not supported as a single number:** the b-value. With mb only it rises steadily with the
  cutoff (1.2–1.3 at mb 4.2, 2.1–2.2 at mb 4.6, `results/robustness/bvalue_vs_cutoff.csv`), so the magnitude
  distribution is not a straight Gutenberg–Richter line in mb (mb saturation is a likely reason).
  Quote b only together with its cutoff, magnitude type and bootstrap interval
  (`results/robustness/bvalues.csv`). Intermediate-depth b-values rest on 3–15 events.
* **Not found:** migrating sequences (0 in every case).

## References

Michailos et al. (2021) Front. Earth Sci.; Him_ides code (GitHub, kemichai); USGS ComCat; ISC Bulletin
(International Seismological Centre); ObsPy; Aki (1965); Utsu (1965); Wiemer & Wyss (2000);
Woessner & Wiemer (2005); Ester et al. (1996, DBSCAN); Chamberlain et al. (2018, EQcorrscan).
