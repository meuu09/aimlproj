# CLAUDE.md: full project context (Him_ides / Himalayan intermediate-depth seismicity)

Read this first. It is the complete handoff, written so a fresh Claude session can continue without redoing work or asking the user to re-explain. Last updated: 2026-11-22.

## 1. The assignment
- University group assignment at IIT Roorkee (group of 7 students). Our topic is #4: **"Analysis of the temporal evolution of Intermediate-depth seismicity beneath the central Himalayas."**
- Professor's requirements:
  - use the source code (Him_ides) and modify it a little;
  - reproduce its results;
  - apply it to a **smaller area, mostly in the Indian Himalaya, using open-source data**;
  - submit online with: Abstract, Methods of data analysis, Results and Discussion, Conclusions, references (including the code).
- **Unresolved:** whether the group submits one report or each student submits their own. An email to the professor was drafted in the group chat; we don't know if it was sent or answered. Ask the user before writing the final report.
- The user is a student, not a seismologist. Explain terms plainly, and ask before choosing anything that changes the science.

## 2. Repository layout
| Path | What it is |
|---|---|
| `Him_ides/` | Vendored copy of https://github.com/kemichai/Him_ides at commit 33acd38 (its `.git` was deleted on purpose so it uploads as ordinary files). It has no LICENSE file; credit the original author (Michailos et al.) in the README and report. Do not edit it; copy and adapt instead. |
| `Him_ides/python/` | `Figures_5-9.py` (main), `temporal_plots.py`, `migration_plots.py`, `Figure_10.py`, `supporting_material/`, `quakeml/` (bundled catalogs), `metadata/` (132 StationXML station files) |
| `Him_ides/gmt/` | GMT shell scripts for map figures (`Figure_1.sh`, `Figure_2.sh`, `Figure_3AB.sh`, `Figure_3CD.sh`, `gmt5.yml`, `him.cpt`, `files/` with station lists, fault `.gmt` files, `.dat` location files) |
| `Him_ides/.ua/` | Understand-Anything knowledge graph (see section 9). Optional. |
| `data/usgs_himachal_kashmir_M3plus_1970-2026.csv` | Our USGS ComCat download (section 4) |
| `him_ides_himachal.py` | **Our adapted analysis script.** Run: `.venv/bin/python him_ides_himachal.py` |
| `results/` | Six figures, `catalog_M4.4_1990plus.csv`, `per_period_summary.csv` |
| `.venv/` | Not in git. Recreate: `python -m venv .venv && .venv/bin/pip install obspy pandas matplotlib scipy numpy` (we used obspy 1.5.2, matplotlib 3.10.9). |

## 3. What Him_ides is
Plotting code for Michailos et al. 2021 (Frontiers in Earth Science) on intermediate-depth (40-150 km) earthquakes beneath the central Himalaya.
- `python/Figures_5-9.py` reads a QuakeML catalog with ObsPy `read_events` from `python/quakeml/`. Catalog read is at lines 33-37.
- Lines 40-44 define the `coords` polygons: NE of Everest 1, SE of one, Udayapur, NW of map, ALL. Lines 63-67 cut the catalog with `tp.cut_cat`.
- `temporal_plots.py` has `get_cat_details`, `cut_cat`, `create_lines`, `interevent_times_single_cluster`, `mag_vs_time_plus_cumulative`, `calc_bvalues`, `multi_mfd_plot`.
- `migration_plots.py` and `Figure_10.py` make the migration and later figures.
- Bundled catalogs: `Diehl_et_al_2017.xml` (919 events) and `Himalayan_Intermediate_depth_earthquake_catalog_Michailos_et_al_2021.xml` (414 events).
- Quirks:
  - `work_dir = os.getcwd()` at line 19, so run from the right folder;
  - hard-coded 2001-2003 dates and Tibet map extents in `temporal_plots.py`;
  - imports `eqcorrscan` (for max-curvature Mc), which is a heavy install.

## 4. Data we got
- **USGS ComCat** (earthquake.usgs.gov FDSN API, count and query endpoints, CSV). Box 31-35°N, 74-79.5°E, 1970-03-05 to 2026-11-14, **566 events, M3.1-6.8** (mostly mb).
- Columns: time, latitude, longitude, depth, mag, magType, nst, gap, dmin, rms, net, id, updated, place, type, horizontalError, depthError, magError, magNst, status, locationSource, magSource.
- We compared four candidate regions at M3, M4 and M5 thresholds and picked the most complete one (Himachal/Kashmir). The catalog is only reasonably complete above about M4.4 (the "M4 completeness" question): below that, the rate rises over time just because more small events are detected.
- Depth caveat: USGS often reports a fixed depth of 10 or 33 km when depth is poorly constrained. See weak point 1.

## 5. Our analysis (`him_ides_himachal.py`)
Parameters at the top:
- `START_YEAR = 1990`
- `MC_CUT = 4.4` (magnitude of completeness cut)
- `DEEP_MIN, DEEP_MAX = 40, 150` (km)
- `FIXED_DEPTHS = (10.0, 33.0)`
- `BOX = dict(lat=(31.0, 35.0), lon=(74.0, 79.5))`

Functions:
- `locator`
- `max_curv(mags, step=0.1)`: re-implements eqcorrscan `calc_max_curv`
- `b_lsq`: the repo's method, a line fit to log10 cumulative counts
- `b_mle`: Aki-Utsu maximum likelihood with binning correction
- `interevent`: inter-event times and coefficient of variation (CV 1 = Poisson, >1 = clustered)

Depth classes: "fixed depth (10/33 km, not measured)", "intermediate depth (40-150 km, measured)", "shallow/other (measured)".

Migration: SVD/PCA of epicentres gives a profile axis (azimuth about 120°); `along_km` is distance along it.

Outputs in `results/`:
1. `Fig1_magnitude_vs_time_cumulative.png`
2. `Fig2_interevent_times.png`
3. `Fig3_magnitude_frequency_bvalue.png`
4. `Fig4_migration_along_arc.png`
5. `Fig5_3D_hypocentres.png`
6. `Fig6_temporal_evolution_summary.png`
7. `catalog_M4.4_1990plus.csv`: 176 events with `depth_class` and `along_km`
8. `per_period_summary.csv`

## 6. Decisions
Ten decisions; six were the user's, four were mine.
- User: target folder `Him_ides/`; install pnpm globally via npm (tooling only); add exclusions `*.xml`, `*.dat`, `*.pyc` to the graph ignore file; **magnitude cut M >= 4.4**; **period 1990 to now**; **make all figures including 3D**.
- Mine (user asked for an explanation and left it to me): region = Himachal/Kashmir box (most complete of the four); flag depths of exactly 10 and 33 km as fixed rather than measured; use both LSQ and MLE b-values; migration along a PCA axis.
- On depth, the user said: "i dontttt think those depths are default values.. if its usgs then we can mark it as reliable." I explained the evidence but did not settle it (weak point 1).

## 7. Results so far (read the caveats in section 8 first)
- 176 events with M >= 4.4 from 1990 onward in the box; depth classes are in `catalog_M4.4_1990plus.csv`.
- b-value: **LSQ 1.82 vs MLE 1.40**. Both are above the usual value of about 1, and the gap is not explained.
- Inter-event time CV and per-period rates are in `per_period_summary.csv` and Fig2/Fig6. Do not claim a real rate change; see weak point 2.

## 8. Weak points: do not overclaim
1. **Fixed-depth interpretation is unverified.** Evidence: 164 events at exactly 33 km with no depth error (none after 2009), and 130 at 10 km with depthError of about 1.9. I did not verify against USGS documentation. The user disagrees and thinks USGS depths can be treated as reliable. Check the docs and resolve before the report states either way.
2. **Rate rise is likely completeness, not real.** More small events are detected over time. Do not interpret the temporal rate increase as physical.
3. **b-value discrepancy uninvestigated** (LSQ 1.82 vs MLE 1.40). LSQ on cumulative counts is known to be biased; the MLE is the safer figure.
4. **Original repo figures not yet reproduced** with the original scripts. This is the professor's first requirement.
5. The Understand-Anything graph lacks the bundled catalogs (the `*.xml` exclusion hid them). Earlier I wrongly said the catalog was only on Zenodo; it is bundled in `python/quakeml/`.
6. Section line ranges for large files in the graph (`temporal_plots.py`, `Figure_10.py`, `plot_catalog_details.py`) are approximate.

## 9. Understand-Anything graph (optional tooling)
- Plugin `understand-anything` was installed from a local tarball marketplace (GitHub clone was too slow). Run via `/understand-anything:understand`.
- `Him_ides/.ua/` holds `knowledge-graph.json` (55 nodes, 85 edges, 5 layers, 12 tour steps), `code-lessons.json` (26 detailed, 12 overview-only), `fingerprints.json`, `meta.json`, `.understandignore`.
- To include the catalogs, add `!python/quakeml/*.xml` to `.understandignore` and rebuild. Not required for the assignment.

## 10. Next steps, in order
1. Ask the user: one group report or individual reports?
2. **Baseline reproduction:** run the original `Him_ides/python/Figures_5-9.py` on the bundled catalog (install `eqcorrscan` or stub its import; run from the right working directory). Save outputs in a separate folder, for example `results_baseline/`.
3. Resolve the depth-flag question and the b-value gap.
4. Optionally repeat for a second region (for example Arunachal) or clip to India's border (the box includes some non-Indian territory).
5. Write the report: Abstract, Methods of data analysis, Results and Discussion, Conclusions, references. Include Michailos et al. 2021, the Him_ides repo, USGS ComCat, ObsPy, Aki 1965 / Utsu for the b-value, Wiemer & Wyss 2000 for Mc, plus the code.

## 11. Working notes
- Use the dedicated tools, keep temporary files out of the project, and don't commit `.venv/`.
- The project is not a git repo until the user uploads it; the user plans to upload to a teammate's private GitHub repo. Never push without being asked.
