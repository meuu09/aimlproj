# Test report

## Commands and results (2026-11-24, branch `audit/final-scientific-validation`)

| Run | Command | Result |
|---|---|---|
| Before the audit | `.venv/bin/python -m pytest -q tests` | 38 passed, 0 failed, 0 skipped |
| After the audit (project venv) | `.venv/bin/python -m pytest -v -W default tests` | **67 passed, 2 xfailed, 0 failed, 0 skipped, 0 warnings** |
| Clean environment | `uv venv --python 3.12 …; uv pip install -r requirements.txt; python -m pytest -q -rs -W default tests` | 67 passed, 2 xfailed, 1 ResourceWarning (unclosed ISC raw file) → fixed in `fetch_isc_bulletin.parse`, rerun: 0 warnings |
| Pipeline, as part of `run_all.py` | step `-m pytest -q tests` | ok |

Total collected: 69 (45 in `test_seis_utils.py`, 8 in `test_data_consistency.py`, 5 functions / 5 tests in
`test_baseline.py`; parametrised cases counted individually).

## Expected failures (strict xfail, not weakened)

| Test | Reason |
|---|---|
| `test_baseline.py::test_published_cluster1_mc` | The published Cluster 1 Mc of 2.85 is not produced by the released code (its 0.5-magnitude bins only allow values ending in .4 or .9) |
| `test_baseline.py::test_published_cluster1_count` | Published Cluster 1 has 248 events; the released polygons select 254 |

`strict=True` means these tests will fail loudly if the discrepancy ever disappears, so the expectation is never silently out of date.

## Test changes made on purpose (not weakening)

`test_depth_class` parameters changed because the classification rule changed (Task 16), from "nominal depth ≤ 150" to "depth + error ≤ 150":

- (150, 5): intermediate → uncertain
- (151, 5): deep → uncertain
- (210, no error): deep → uncertain

(145, 5) → intermediate and (156, 5) → deep were added to pin the new boundaries.

## What the tests cover

| Area | Tests |
|---|---|
| Completeness | `max_curv` binning, eqcorrscan Mc identical to the verbatim 0.4.2 code, GFT and MBS recover synthetic Mc (synthetic fixtures are labelled as such), narrow-window NaN |
| Gutenberg-Richter | exact G-R LSQ fixture, MLE formula and error, MLE recovers b = 1 on 20 000 synthetic events, Shi-Bolt formula, bootstrap repeatability, reliability labels, constant-magnitude limit log10(e)/(Δm/2) |
| Depth | interval classification at both bounds, missing errors, NaN depth |
| Rates and time | Garwood CIs against tabulated values (n = 0, n = 10), CV test (regular → low, burst → high), dispersion index, inter-event ordering and units |
| Distances and geometry | haversine 1° = 111.19 km, 3-D Pythagoras, PCA azimuth of a known 130° line, DBSCAN [lat, lon] order |
| Catalog matching | one-to-one, closest first, competing candidates, transitive USGS-ISC-USGS chain, ISC-only groups, pair-level flags, far pairs |
| Parsers and data | ISF parser field by field on two genuine 2005 ISC records; whole-year parser regression; USGS depth metadata belongs to the CSV origin; combined rows use one whole origin; no sentinel errors; dates and box; every USGS event kept once |
| Migration | all candidates kept, overlap rule, inclusive window boundary, linear and diffusive fits recovered exactly, permutation p |
| Monte Carlo | reproducible with a fixed seed; missing depth errors are not perturbed |
| Baseline | released-code values; published values as strict xfail; Cluster 1 b at the published Mc within 1σ |
| Edge cases | empty, single, constant and NaN inputs return NaN or empty results, not crashes |

## Tests that need data

Every data test reads files that are committed to the repository: `data/`, plus `results_baseline/` for `test_baseline.py`, which skips if that folder is missing. No test needs a network connection. CI (`.github/workflows/tests.yml`) runs all of them plus `run_all.py --quick`. It has not yet run on GitHub because the branch has not been pushed.
