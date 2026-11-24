# Initial audit (Stage A)

Branch `audit/final-scientific-validation`, created from `audit/qc-and-fixes` at commit e2e7825, on 2026-11-24.
This stage changed no code, data or results. The per-task evidence is in `TASK_STATUS.csv` (85 rows).

## 1. Baseline established before any change

| Check | Result |
|---|---|
| Test suite | `.venv/bin/python -m pytest -q tests`: **38 passed**, 0 failed, 0 skipped |
| Full pipeline | `.venv/bin/python run_all.py --baseline`: every step ok |
| Reproducibility | A rerun gives byte-identical files in `results/` and `results_baseline/`, figures included (git shows no change) |
| Upstream | `Him_ides/` matches kemichai/Him_ides at 33acd38. The GitHub API reports `license: null`. |

## 2. Status summary

| Status | Tasks |
|---|---|
| COMPLETE | 7 (5, 21, 33, 47, 80, 81, 82) |
| PARTIALLY COMPLETE | 44 |
| INCORRECT | 10 (6, 8, 11, 12, 23, 52, 55, 65, 73, 74) |
| NOT IMPLEMENTED | 23 |
| BLOCKED | 1 (53: Omori; there are no aftershock sequences to fit) |
| NOT APPLICABLE | 0 |

The waveform methods (EQTransformer, template matching, kpick, NonLinLoc, HypoDD) are out of scope, as the brief says, so no task depends on them.

## 3. Source caveat: the paper

The paper page was fetched, but only the abstract, introduction and Sections 2.1-2.4 came back as verbatim text. The published values used below are therefore **not yet verified from the paper's own text**: Cluster 1 = 248 events, Mc 2.85, b 0.96, the 20-day / 20-km 3-D migration window, and the figure captions. The parts confirmed verbatim are: Mc 2.4 (abstract), 414 events, ">40 km" for intermediate depth, observed depths of 40-110 km, and the dehydration interpretation. **Tasks 1-3 and 56 need the PDF.**

## 4. Confirmed findings (evidence reproduced in this stage)

### Original-paper reproduction
- **Task 1 (Cluster 1 Mc).** The released `calc_max_curv` uses 0.5-magnitude bins, so Mc + 0.4 can only end in .4 or .9. On the bundled catalog, Cluster 1 (254 events) gives Mc = 2.0 + 0.4 = **2.4**. This is the same whether events are selected by the last, first or preferred origin. Bin widths of 0.2, 0.1 and 0.05 give 3.4, 2.9 and 3.25. **No setting of the released algorithm gives 2.85.** The value "2.43" in our earlier notes is the smallest magnitude used in the fit (the figure legend label), not Mc. This is a reporting error on our side; the calculation itself is correct.
- **Task 2 (b-values).** The released least-squares fit gives Cluster 1 b = 0.897 at the released Mc. A maximum-likelihood fit at M ≥ 2.85 gives 0.98 ± 0.13 (N = 55). The gap comes from the threshold and the event selection (254 vs 248), not from our wrapper, which calls the original function. The authors' catalog also mixes magnitude labels (ML 397, ML* 10, M_init 7); leaving out M_init does not change Mc.

### Catalog construction
- **Task 8 (critical).** 187 combined events pair a USGS epicentre with an ISC-agency depth. ISC-only events can take the prime epicentre from one agency and the depth from another. A hypocentre is a joint solution, so these rows are not internally consistent.
- **Task 11.** In duplicate grouping, a chain USGS A, then ISC X, then USGS B can merge two independent USGS events and flag one of them. This is latent: today all 7 flagged rows are ISC-only, so no data are affected yet.
- **Task 12.** The `flagged` column in `close_pairs.csv` uses event-level status, not whether that pair meets the rule. Also latent: today 0 of 11 pairs are misreported.
- **Task 15.** Equivalent magnitude labels are not normalised (MW 261, mww 10, mwc 3, mw 4 ...; ML 906 vs mL 10). The "mb only" filters do match mb exactly.
- **Task 5 (complete).** For all 566 USGS events the depth-type metadata matches the CSV origin (depth within 0.5 km, same location source).

### Depth
- **Task 16.** The intermediate rule applies the uncertainty at 40 km but not at 150 km. Three combined events are "intermediate" even though their error bars cross 150 km.
- **Tasks 18 and 39.** The Monte Carlo treats a missing depth error as zero (`fillna(0)`), i.e. as a perfect depth. The depth-class rule itself handles missing errors correctly (they become "uncertain").
- **Task 19.** None of the solved depths has a negative or zero error, but 137 combined solved events have depth errors above 50 km. These have not been examined.

### Completeness and b-values
- **Task 23 (critical).** Verified: with the authors' method, USGS Mc is **4.9 in the 1970s, 1980s and 1990s**, 3.9 in the 2000s and 4.4 in the 2010s and 2020s. The per-period rates use M ≥ 4.4 from 1990 onward, so **the 1990s rate comes from an incomplete catalog**. M 4.4 should be presented as an exploratory threshold.
- **Task 30.** The main figures and CSVs show MLE b-values for groups of 5-15 events without any reliability label.
- **Task 31.** CLAUDE.md attributes the mb b-value rise to "mb saturation likely". That cause has not been tested.

### Temporal analysis and migration
- **Task 43 (critical, research objective).** There is no dedicated temporal summary for the intermediate-depth events. `per_period_summary.csv` covers all events only.
- **Task 52.** The Båth's-law difference is computed between the two largest intermediate events over 36 years, which are unrelated. That is not a valid use of the law.
- **Task 55.** Candidate sequences whose events are a subset of an earlier sequence are thrown away before saving.
- **Task 53.** Blocked: there are 0 credible sequences, so an Omori fit is not possible.

### Visualisation and pipeline
- **Task 65.** The 3-D plot has degrees on the horizontal axes and km on the depth axis, which distorts the geometry.
- **Task 73.** `run_all.py` does not stop on failure. If the catalog merge fails, the combined analyses still run on the old CSV.
- **Task 74.** `--quick` is documented as "USGS only", but it still runs the steps that need the combined catalog.
- **Task 6.** Failed USGS metadata requests are never retried. No impact today: there are 0 failed rows.

### Licensing
- **Task 85 (critical, legal).** Upstream has no licence, which by default means all rights reserved. **`meuu09/aimlproj` is public** and contains a copy of `Him_ides/` (code and catalogs). Attribution is present, but attribution alone does not grant the right to redistribute.

## 5. Decisions needed from you before Stage B

1. **Origin consistency (Task 8).** Option (a): use one whole origin per event, chosen by agency priority (ISC-EHB > ISC > NEIC > NDI > USGS), and keep the other origins in a side table. Option (b): keep the current mixing and document it. I recommend (a). It will change the intermediate-event counts.
2. **Thresholds (Tasks 23 and 28).** Keep M 4.4 as the exploratory reference and add M 4.9 (conservative, whole period) and post-2000 runs, then choose the headline from the completeness evidence. I recommend this.
3. **Depth rule (Task 16).** Apply the error at both bounds (depth − error ≥ 40 **and** depth + error ≤ 150). I recommend this; it moves 3 combined events to "uncertain".
4. **Study area (Task 14).** Keep the box, documented as a geological region that extends beyond Indian-administered territory, or add an India-only sensitivity run.
5. **Licensing (Task 85).** Ask the author for permission, or remove `Him_ides/` from the public repository and have `reproduce_baseline.py` clone it at commit 33acd38. I recommend asking the author and making the repository private until then. This affects your teammate's GitHub repository, so it is your call.
6. **Paper PDF.** Tasks 1-3 and 56 need the paper's exact text, which the web fetch did not return. If you can save the PDF into the project, I will verify the values from it.

## 6. Proposed Stage B order (after approval)

1. Baseline: comparison table, correct the Mc labelling, baseline tests (Tasks 1-4, 77).
2. Catalog: origin consistency, transitive duplicate fix, pair-level QC, magnitude-type normalisation, dates, retries (Tasks 6, 8, 11-13, 15, 75).
3. Depth: symmetric rule, missing-error handling, error audit, 40/50/60/70 km sensitivity (Tasks 16-20).
4. Completeness: threshold comparison, Mc per catalog, type and cluster, sample-size labels (Tasks 23-31).
5. Temporal: intermediate summary, Poisson confidence intervals, CV against simulated Poisson, Båth only within sequences (Tasks 43-52).
6. Clustering: missing-error handling in the Monte Carlo, joint propagation, 3-D DBSCAN, k-distance plot, polygon tests (Tasks 33-42).
7. Migration: keep all candidates, region-wide search, window sensitivity (Tasks 54-63).
8. Visualisation: 3-D plot in km, cross-sections (Tasks 64-66).
9. Pipeline: fail-fast, `--quick`, run manifest, CI (Tasks 72-74, 78-84).
10. Report documents (Tasks 17, 67-71) and the final deliverables.

Each step will be validated with tests, a rerun, and a before/after comparison (`compare_results.py`).
