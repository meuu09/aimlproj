# Temporal evolution of intermediate-depth seismicity beneath the NW (Indian) Himalaya

IIT Roorkee group assignment (topic 4). We reproduce the catalog-based analysis of **Michailos, Carpenter & Hetényi (2021)**, *Spatio-Temporal Evolution of Intermediate-Depth Seismicity Beneath the Himalayas*, Front. Earth Sci. 9:742700 (doi:10.3389/feart.2021.742700). We use their code ([kemichai/Him_ides](https://github.com/kemichai/Him_ides)), adapt it, and apply it to the Himachal–Kashmir Himalaya (31–35°N, 74–79.5°E) with open catalogs from **USGS ComCat** and the **ISC Bulletin**.

The full audit of methods, code and results is in [`audit/FINAL_AUDIT.md`](audit/FINAL_AUDIT.md).

## Research objective

Do intermediate-depth earthquakes (here 40–150 km, with their depth error) beneath the NW Himalaya cluster in space, change in rate, occur in sequences, or migrate, as reported for the central Himalaya? And how much of any answer is limited by the catalogs themselves?

## Setup

Python 3.12 is required. The original code imports eqcorrscan, which does not install on 3.12, so `reproduce_baseline.py` supplies the one function it needs (`calc_max_curv` from eqcorrscan 0.4.2).

```bash
uv venv --python 3.12 .venv          # or: python3.12 -m venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python -m pytest -q tests
.venv/bin/python run_all.py --baseline
```

**Windows (PowerShell):** `py -3.12 -m venv .venv`, then `.venv\Scripts\python -m pip install -r requirements.txt` and the same commands with `.venv\Scripts\python`.

`run_all.py` rebuilds every result from `data/` in dependency order (about a minute, no downloads). It skips every step that depends on a failed step, logs to `results/run_log.txt`, runs the tests, and writes `results/MANIFEST.json` (checksums, versions, seeds, commit).

Options:
- `--baseline`: also rerun the original code first;
- `--quick`: USGS-only analyses;
- `--download`: refresh the downloads.

Two runs give identical outputs (see `audit/REPRODUCIBILITY.md`).

## Original methodology vs ours

Details are in [`audit/PAPER_METHOD_COMPARISON.md`](audit/PAPER_METHOD_COMPARISON.md).

- **The paper:** waveforms from temporary networks (2001–2003), ML magnitudes and HypoDD relocations, 414 events. Mc by maximum curvature + 0.4. b-value from a least-squares fit to the cumulative distribution. Inter-event times. A migration search within 20 days and 20 km (3-D) of any first event.
- **This project:** bulletin catalogs only (no waveform processing), 1970–2026, analysis from 1990. The paper's Mc, b-value, inter-event and migration methods are reproduced. We add:
  - depth classes that use the depth error;
  - completeness by catalog, decade, magnitude type and cluster with five methods;
  - maximum-likelihood b-values with sample-size rules;
  - exact Poisson rate intervals and Poisson tests;
  - DBSCAN (2-D and 3-D) with location-error Monte Carlo;
  - a region-wide migration search with window sensitivity;
  - Global CMT focal mechanisms.

## Data and quality control

| Data | Script | Notes |
|---|---|---|
| USGS ComCat events (M ≥ 3) | – | `data/usgs_himachal_kashmir_M3plus_1970-2026.csv` |
| USGS depth-type flag (solved vs fixed) | `fetch_depth_type.py` | 566/566 checked against the CSV origin; failed requests are retried |
| ISC Bulletin (ISF, all agencies) | `fetch_isc_bulletin.py` | One **whole origin** per event: the first trusted agency (ISC-EHB > ISC > NEIC > NDI) with a free depth and an error, else the ISC prime origin |
| USGS + ISC combined | `build_combined_catalog.py` | Matched on the ISC prime origin (20 s, 100 km, one-to-one). Each row is one self-consistent origin; replaced USGS origins are kept in `usgs_*` columns. `magType_norm` groups labels of the same scale. Possible duplicates are flagged (ISC-only rows only) |
| Global CMT solutions | `fetch_focal_mechanisms.py` | 24 events; P/T/B axes and World Stress Map regime |
| GEM Global Active Faults subset | – | `data/gem_active_faults_box.geojson`, CC BY-SA 4.0 (Styron & Pagani 2020) |
| Quality reports | `catalog_qc.py` | `results/qc/QC_REPORT.md`: filter funnel, magnitude types, depth quality by decade, uncertainties by agency, match-threshold sensitivity, review of matches > 50 km, duplicate pairs |

**Study area.** The box is a geological study region. It includes areas administered by Pakistan and China, and is not clipped to political borders.

## Definitions

* **Depth classes** (from the whole error bar):
  - intermediate: 40 ≤ depth − error and depth + error ≤ 150 km;
  - shallow: depth + error < 40 km;
  - deep: depth − error > 150 km;
  - uncertain: otherwise, or no error (a missing error is never treated as zero).
* **Solved depth:** USGS depth-type "from location" or "constrained by depth phases"; for ISC, a free depth with an error from a trusted agency. Two cases are run throughout:
  - `excluding_fixed` (main case);
  - `including_fixed`.
* **Magnitudes are never converted** between scales. mb-only analyses give a homogeneous scale.
* **Mc:** the authors' method (eqcorrscan 0.4.2 maximum curvature + 0.4), plus the Wiemer–Wyss mode, GFT and MBS.
* **M 4.4 is an exploratory threshold.** USGS is complete only at M 4.9 in the 1990s.
* **b-value:** Aki–Utsu maximum likelihood (with the paper's least-squares fit for comparison). Not reported for N < 20; "rough" for 20 ≤ N < 50.

## Scripts and outputs

| Script | Output |
|---|---|
| `reproduce_baseline.py`, `baseline_comparison.py` | `results_baseline/` (Figures 5–10, `baseline_comparison.md`) |
| `him_ides_himachal.py`, `kishtwar_cluster.py`, `depth_clusters.py` | `results/<catalog>/<case>/` (options `--include-fixed`, `--catalog combined`, `--drop-duplicates`) |
| `temporal_summary.py` | `results/temporal/`: rates with CIs, CV tests and b per depth population and completeness scenario; `intermediate_summary.csv` |
| `migration_search.py` | `results/migration/`: all candidate sequences, window sensitivity |
| `clustering_checks.py` | `results/clustering/`: 2-D vs 3-D DBSCAN, k-distance, polygon variants, PCA strikes, space-time feasibility, `base_map.png` |
| `mc_through_time.py`, `bvalue_robustness.py`, `cluster_robustness.py`, `depth_sensitivity.py` | `results/robustness/` |
| `catalog_qc.py` | `results/qc/` |

Shared functions are in `seis_utils.py`; the tests (67, plus 2 documented expected failures) are in `tests/`.

## Verified results

1. **Baseline reproduction** (`results_baseline/baseline_comparison.md`).
   - Full catalog: N = 414, Mc 2.4 and b 0.917 (published 0.91) are reproduced with the released code.
   - Cluster 1 b matches at the published Mc: 0.942 (least squares) and 0.976 ± 0.13 (maximum likelihood) against 0.96.
   - Not reproduced: the published Cluster 1 Mc of 2.85 cannot come from the released code (0.5 bins), and Cluster 1 has 254 events against the published 248. The paper's migration-search code is not public.
2. **Intermediate-depth events are few.** Solved depths, since 1990: 25 (USGS) and 50 (combined) events; at M ≥ 4.4, 8 and 15.
3. **A concentration beneath Kishtwar–Chamba** is found by the hand-drawn outline and by DBSCAN in 2-D and 3-D (adjusted Rand index 0.93–1.0). It is **fragile**:
   - it shrinks from 17 to 5 events (combined) when the lower depth bound moves from 40 to 50 km;
   - its membership changes with small outline changes;
   - in the Monte Carlo, at least 90 % of its events stay in one cluster in only 68–79 % of runs, and in 4–15 % if the errors are twice as large;
   - the best-matching cluster contains many unrelated events.
4. **Depth reliability is the main limitation.**
   - The regional crust is about 44–54 km thick (Sharma et al. 2018).
   - All 24 Global CMT centroids are at 41 km or shallower.
   - The one catalog "intermediate" event with a GCMT solution (49.9 km) has a centroid at 25 km.
   - So these events cannot be assigned to the lower crust or the mantle, and some may be mislocated in depth.
5. **Rates and clustering in time.** Intermediate-event inter-event times are not distinguishable from a stationary Poisson process under any completeness scenario. Rate changes between decades are not interpretable: completeness and the share of solved depths both change.
6. **No sequences and no migration at intermediate depth** under any window from 5 days/10 km to 60 days/50 km. This reflects detection limits (M ~4, location errors comparable to the window); it does not show that migration is absent.
7. **b-values depend on the magnitude cutoff and the magnitude type.** The distributions do not follow a clean Gutenberg–Richter law (the GFT and MBS methods mostly find no Mc). The cause is not determined.

## Scientific limitations

The list is in `audit/FINAL_AUDIT.md` §3 and `audit/SCIENTIFIC_ASSUMPTIONS.md`. In short:
- bulletin depths with errors of about 8–15 km;
- mixed magnitude scales;
- completeness changing over time;
- small samples;
- uncertainty definitions not verified;
- no relocation.

Metamorphic dehydration and eclogitization (the paper's interpretation) are hypotheses that this catalog cannot test.

## Licensing and attribution

- **Him_ides** (Michailos et al.) has **no licence** file. We do not redistribute it: `reproduce_baseline.py` clones it at commit 33acd38 when needed, and the original authors are credited here and in the report. Earlier commits of this repository contained a copy, and permission from the author has not been obtained.
- **GEM faults:** CC BY-SA 4.0 (attribution inside the file).
- **Global CMT:** Dziewonski et al. (1981); Ekström et al. (2012).
- **ISC Bulletin:** International Seismological Centre (cite per isc.ac.uk/citations).
- **USGS ComCat:** public domain.

## References

Michailos, Carpenter & Hetényi (2021) Front. Earth Sci. 9:742700; Him_ides (GitHub, kemichai); USGS ComCat; ISC Bulletin; Global CMT (Dziewonski et al. 1981; Ekström et al. 2012); Styron & Pagani (2020, GEM faults); ObsPy; Aki (1965); Utsu (1965); Shi & Bolt (1982); Wiemer & Wyss (2000); Woessner & Wiemer (2005); Cao & Gao (2002); Garwood (1936); Båth (1965); Shapiro et al. (1997); Ester et al. (1996, DBSCAN); Zoback (1992); Chamberlain et al. (2018, EQcorrscan); Sharma et al. (2018, EGU2018-11932); Shamim et al. (2024, EGU24-15623). Full list: `audit/SCIENTIFIC_ASSUMPTIONS.md`.
