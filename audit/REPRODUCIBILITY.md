# Reproducibility

## Environment

Python 3.12 (verified with 3.12.16 in a fresh environment). Exact versions are pinned in `requirements.txt`:
numpy 2.5.5, pandas 3.1.1, matplotlib 3.11.2, scipy 1.18.1, scikit-learn 1.9.1, shapely 2.2.0, obspy 1.5.2 (baseline only), pytest 9.2.0. `pyproj` was removed because no script uses it.

```bash
uv venv --python 3.12 .venv                       # or: python3.12 -m venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python -m pytest -q tests
.venv/bin/python run_all.py --baseline            # full rebuild, ~1 min
```

Windows: use `.venv\Scripts\python` with the same arguments.

## Execution modes (`run_all.py`)

| Mode | Steps |
|---|---|
| default | ISC parse (offline) → combined catalog → QC → completeness → focal-mechanism parse → 18 analysis runs (3 catalogs × 2 depth cases × 3 scripts) → temporal → migration → clustering → b-value, cluster and depth robustness → tests → manifest |
| `--baseline` | adds, first: the original Him_ides code (cloned at commit 33acd38 if `Him_ides/` is missing) and the published-value comparison |
| `--quick` | USGS-only analyses (6 runs) → tests → manifest; no combined catalog, works offline |
| `--download` | also refreshes the USGS depth types, the ISC bulletin (missing years only) and Global CMT |

If a step fails, every step that depends on it is skipped and logged in `results/run_log.txt`, and the run exits with code 1. This was verified by simulating a failure of `build_combined_catalog.py`.

## Data provenance

| File | Source | Retrieved |
|---|---|---|
| `data/usgs_himachal_kashmir_M3plus_1970-2026.csv` | USGS ComCat FDSN event service, box 31-35°N 74-79.5°E | 2026-11-23 |
| `data/usgs_depth_type.csv` | ComCat per-event GeoJSON, preferred origin `depth-type` (566/566 match the CSV origin) | 2026-11-23 |
| `data/isc_raw/isc_isf_*.txt` | ISC Bulletin, COMPREHENSIVE, ISF (IMS1.0), per year 1970-2026 | 2026-11-23 |
| `data/gcmt_box_cmtsolution.txt` | Global CMT web search (CMTSOLUTION) | 2026-11-24 |
| `data/gem_active_faults_box.geojson` | Subset of GEM Global Active Faults (CC BY-SA 4.0; attribution in the file) | 2026-11-24 |
| Upstream code | github.com/kemichai/Him_ides @ 33acd38 (no licence; not redistributed) | cloned on demand |

The SHA-256 checksum of every input file is in `results/MANIFEST.json`.

## Random seeds

All random numbers come from `numpy.random.default_rng` with fixed seeds:

| Function | Seed | Draws |
|---|---|---|
| b-value bootstrap | 0 | 2000 or 1000 resamples |
| Mc sliding-window bootstrap | 0 | 200 resamples |
| CV Poisson test | 0 | 2000 simulations |
| Migration permutation test | 0 | 2000 permutations |
| Cluster Monte Carlo | 0 | 500 runs per catalog and scenario |

## Reproducibility test

`run_all.py --baseline` was run twice on the same commit (2026-11-24) and the two `MANIFEST.json` files compared:
- all **93 output tables (CSV/MD) have identical SHA-256** checksums;
- `git status` showed no figure (PNG/JPG) changed;
- the input checksums are identical.

Only `MANIFEST.json` (date) and `run_log.txt` (timings) differ, as expected.

Two earlier runs before the audit were also byte-identical (66 CSVs).
