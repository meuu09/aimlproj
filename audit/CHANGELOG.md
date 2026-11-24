# Change log (branch `audit/final-scientific-validation`)

Validation for every entry: the full test suite (67 passed, 2 xfailed) and `run_all.py --baseline` (33 steps ok). Numerical changes are listed in `BEFORE_AFTER_RESULTS.csv` (443 rows, 190 changed).

| File | What changed | Why | Tasks | Numbers changed? |
|---|---|---|---|---|
| `seis_utils.py` | Interval depth rule at both bounds; `CATALOG_START`/`CATALOG_END` applied in `load_catalog`; `norm_magtype`; `flag_duplicates` rewritten (no USGS chains, pair-level flags); new `b_reliability`, `mc_gft`, `mc_mbs`, `poisson_ci`, `cv_poisson_test`, `dispersion_index`, `hypo_dist_km`, `candidate_sequences`, `resolve_overlaps`, `migration_fit` | Confirmed bugs and missing methods | 11-13, 15-16, 27, 30, 46, 49-50, 55-59, 61 | Yes: depth classes (3 combined events), date filter |
| `fetch_isc_bulletin.py` | One whole origin per event (trusted agencies only), prime origin kept for matching; parses time error, rms, smaj, smin, nsta, gap; `--parse-only`; skips years already downloaded; closes files | Mixed hypocentres; dropped uncertainty fields | 7, 8 | Yes: 195 ISC event coordinates |
| `fetch_depth_type.py` | Retries failed rows, back-off, failure log, one row per event | Failures were never retried | 6 | No (0 failures; file byte-identical) |
| `build_combined_catalog.py` | Rewritten: match on ISC prime origin; whole-origin rows; `usgs_*` alternative columns; `location_source`, `origin_id`, `magType_norm`; sentinel errors set to missing | Task 8 | 8, 15, 19 | Yes: combined solved 1095 → 1086 (1990+), intermediate 55 → 50 |
| `catalog_qc.py` | Pair-level duplicate eligibility and flags; new tables: uncertainty by source, USGS location quality, match-threshold sensitivity, review of matches > 50 km | Misleading QC; missing diagnostics | 9, 10, 12, 19, 22 | New tables |
| `reproduce_baseline.py` | Clones upstream at 33acd38 if missing; records the released-code Mc and the magnitudes; quotes the CSV field | Licensing; mislabelled Mc; malformed CSV | 1, 85 | Baseline CSV columns (figures unchanged) |
| `baseline_comparison.py` (new) | Published vs reproduced table (CSV + MD) | Missing deliverable | 1-4 | New |
| `him_ides_himachal.py` | b suppressed below N = 20 and labelled; rate CIs; intermediate count per period; 3-D plot in km with stated ×2 exaggeration; new class labels | Small-sample b; misleading 3-D geometry | 30, 45-46, 65 | per_period b for N < 20 removed |
| `kishtwar_cluster.py` | Shared sequence search (all candidates kept, overlaps marked); linear/diffusive fits with errors; Båth only within sequences (paper-style value labelled); b labels; missing-error markers | Discarded candidates; misused Båth | 51-52, 55-59 | Summary keys renamed |
| `depth_clusters.py` | Class labels for the new rule | Task 16 | 16 | Labels |
| `migration_search.py` (new) | Region-wide search, window grid, resolvability, permutation test | Polygon-only search | 54, 57, 60-63 | New: 0 intermediate candidates |
| `temporal_summary.py` (new) | Completeness scenarios × depth populations × periods: rates and CIs, CV test, dispersion, b if N allows, completeness flag | No intermediate-depth temporal summary | 23, 28, 43-50 | New |
| `mc_through_time.py` | Mc by five methods for every catalog, type and the Kishtwar outline; `magType_norm` | Single-method, single-catalog Mc | 24-27, 29 | New table |
| `bvalue_robustness.py` | Adds the duplicate-removed catalog; `magType_norm` | Task 79 | 79 | New rows |
| `cluster_robustness.py` | Joint depth+epicentre Monte Carlo; missing depth errors not perturbed; error-interpretation scenarios; 50/70/80/90 % recovery, contamination, fragments; duplicate-removed catalog | Missing errors were set to zero; separate, looser tests | 18, 39-41, 79 | Yes: band recovery and Monte Carlo tables |
| `clustering_checks.py` (new) | 3-D DBSCAN, 2-D vs 3-D, k-distance, polygon variants, PCA strikes, space-time feasibility, base map | Missing checks | 33-38, 42, 64 | New |
| `depth_sensitivity.py` (new) | Lower bounds 40-70 km, upper bounds 100/150 km | Task 20 | 20 | New |
| `fetch_focal_mechanisms.py` (new) | Global CMT download, P/T/B axes, WSM regime, catalog match | Task 68 | 68 | New |
| `run_all.py` | Dependency-ordered steps, skip on upstream failure, log, real `--quick`, `--download`, tests and manifest as final steps | No fail-fast; wrong quick mode | 72-74 | — |
| `write_manifest.py` (new) | Checksums, versions, seeds, commit | Task 83 | 83 | — |
| `.github/workflows/tests.yml` (new) | CI: install, tests, quick run | Task 84 | 84 | — |
| `.gitignore` / `Him_ides/` | Upstream copy no longer tracked on this branch | No upstream licence | 85 | — |
| `requirements.txt` | `pyproj` removed (unused) | Task 80 | 80 | — |
| `tests/` | +31 tests (see TEST_REPORT.md) | Tasks 75-78 | 75-78 | — |
| `data/*` | Regenerated ISC event table and combined catalog; new GCMT text and GEM subset; raw ISF unchanged | Tasks 8, 64, 68 | — | Yes (above) |
| `results/`, `results_baseline/` | Regenerated by the corrected pipeline; new `temporal/`, `migration/`, `clustering/`, `geology/` folders | — | — | See BEFORE_AFTER_RESULTS.csv |
| `README.md`, `CLAUDE.md` | Updated to the corrected methods and results | Task deliverable 10 | — | — |
| `audit/*` | Audit deliverables; `collect_metrics.py` and `build_before_after.py` regenerate the before/after table | — | — | — |
