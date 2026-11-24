# Himalayan intermediate-depth seismicity

This repository adapts [Him_ides](./Him_ides/README.md), the plotting code for Michailos et al. (2021), to a Himachal–Kashmir study box. The box is 31–35°N, 74–79.5°E and crosses national boundaries; it is a western-Himalaya study area, not the central Nepal–southern Tibet footprint of the original paper. State that scope explicitly in any report.

## Repository map

| Path | Purpose |
|---|---|
| `data/usgs_himachal_kashmir_M3plus_1970-2026.csv` | Stored USGS ComCat M3+ snapshot. The original query was not previously recorded as a script. |
| `data/usgs_depth_type.csv` | Per-event ComCat depth-type metadata for the stored USGS snapshot. |
| `data/isc_raw/` | Annual ISF responses downloaded from ISC for 1970–2026. The raw files are retained. |
| `data/isc_bulletin_events.csv`, `data/isc_bulletin_origins.csv` | Original parsed ISC event/origin tables. |
| `data/combined_usgs_isc_catalog.csv` | Original combined catalog; retained for comparison. |
| `data/*_v2.csv` | Versioned parsed/combined tables with origin and magnitude provenance. They do not replace the original tables. |
| `fetch_usgs_catalog.py`, `fetch_depth_type.py`, `fetch_isc_bulletin.py`, `build_combined_catalog.py` | Collection, metadata, parsing, and catalog combination. |
| `catalog_qc.py`, `catalog_quality_report.py` | Match diagnostics and non-destructive catalog quality reports. |
| `him_ides_himachal.py` | Regional temporal, magnitude-frequency, PCA-profile, and 3D analyses. |
| `kishtwar_cluster.py` | Kishtwar–Chamba polygon analysis and fixed-window sequence search. |
| `depth_clusters.py` | Spatial DBSCAN clusters separated by depth-quality group. |
| `reproduce_baseline.py` | Runs the original Him_ides Python Figures 5–9 without editing that source tree. |
| `results/` | Existing results plus newly generated quality reports and v2 analyses. |
| `results_baseline/` | Existing reproduction of original Python Figures 5–9 and b-value table. |
| `Him_ides/` | Vendored original source, its bundled QuakeML catalogs, StationXML, GMT scripts, and map inputs. Do not edit this tree. |

The original QuakeML catalogs contain 414 Michailos et al. events and 919 Diehl et al. events. The baseline runner uses the 414-event catalog. It has not run original Figure 10 or the GMT map/cross-section scripts.

## Windows setup

Run these commands in PowerShell from the repository root:

```powershell
py -3.12 -m venv .venv
& .\.venv\Scripts\python.exe -m pip install --upgrade pip
& .\.venv\Scripts\python.exe -m pip install -r requirements.txt
& .\.venv\Scripts\python.exe -m unittest discover -s tests -v
```

`scikit-learn` is required by DBSCAN. ObsPy is needed by the original catalog-detail and Figure 10 tools. The baseline wrapper supplies the original `eqcorrscan 0.4.2` maximum-curvature routine, so installing `eqcorrscan` is not required. The original GMT scripts require GMT classic and Unix shell utilities; their display commands are not a native PowerShell workflow.

## Pipeline

1. `fetch_usgs_catalog.py` records an explicit ComCat query and writes a new CSV plus a JSON query sidecar. It refuses to overwrite an existing file unless `--force` is provided.
2. `fetch_depth_type.py` queries ComCat origin products and writes a versioned depth-flag table. Failed requests remain marked and can be retried on a later run.
3. `fetch_isc_bulletin.py --parse-only` parses the stored annual ISF files without contacting ISC. Its default output is versioned. ISF depth errors retain their 90% confidence label.
4. `build_combined_catalog.py` uses the existing time/epicenter rule (≤20 s and ≤100 km) with a one-to-one greedy assignment. It writes a versioned catalog and a diagnostic row for each accepted match. Review flags do not reject matches.
5. `catalog_quality_report.py` writes source, magnitude-type, depth-quality, filtering, and origin-consistency reports under `results/catalog_quality/`.
6. The analysis scripts use an M3 input floor for comparisons with the stored USGS M3+ snapshot. The original combined catalog remains available through `--catalog combined`; the provenance-rich table uses `--catalog combined-v2`.

Example PowerShell commands:

```powershell
& .\.venv\Scripts\python.exe fetch_isc_bulletin.py --parse-only
& .\.venv\Scripts\python.exe build_combined_catalog.py
& .\.venv\Scripts\python.exe catalog_quality_report.py
& .\.venv\Scripts\python.exe him_ides_himachal.py --catalog combined-v2
& .\.venv\Scripts\python.exe him_ides_himachal.py --catalog combined-v2 --include-fixed
& .\.venv\Scripts\python.exe kishtwar_cluster.py --catalog combined-v2
& .\.venv\Scripts\python.exe depth_clusters.py --catalog combined-v2
& .\.venv\Scripts\python.exe reproduce_baseline.py
```

To run DBSCAN on every event with a reported magnitude, use `--min-mag 0`. The standard M3+ comparison is the default. These choices are recorded in the plot stamp and output tables.

## Scientific interpretation and limits

- The magnitude value is always paired with its original `magType` and source. This project performs no magnitude conversions. An M4.4 cutoff across `mb`, `mL`, `md`, and `MW` is not a homogeneous threshold.
- The ISC `depthError` is a 90% error in kilometres; the ComCat `depthError` is retained with its source and is not assigned the ISC confidence level. Missing errors stay missing in plots and reports. Very large reported values are flagged for review, not silently removed.
- ISC event-level summaries can use prime-origin time/coordinates and a depth from another agency origin. The v2 tables retain both origin identities and coordinates so this can be checked.
- The M3 floor is applied before decade counts and frequency-magnitude summaries. The reports also show how many rows below M3 were present in the combined source.
- Maximum-curvature estimates are shown by magnitude type/source. The original paper's +0.4 correction is retained for reproduction, but estimates from a mixed-scale catalog are not combined into one b-value.
- Period rates, b-values, PCA projections, DBSCAN groups, and Kishtwar sequence searches are descriptive or exploratory until catalog completeness, source homogeneity, location uncertainty, and sample size are assessed. A projected PCA axis is not itself evidence of earthquake migration.
- The current data are not enough to support a physical rate-change, depth-migration, or tectonic-cause claim. A zero count under the chosen 20-day/20-km sequence rule is specific to that rule and catalog selection.

## References

- Michailos, K., Carpenter, N. S., and Hetényi, G. (2021), [Spatio-Temporal Evolution of Intermediate-Depth Seismicity Beneath the Himalayas](https://doi.org/10.3389/feart.2021.742700).
- [Him_ides source repository](https://github.com/kemichai/Him_ides) and its [Python instructions](./Him_ides/python/README.md).
- [USGS ComCat documentation](https://www.usgs.gov/programs/earthquake-hazards/anss-comprehensive-earthquake-catalog-comcat-documentation) and [FDSN Event Web Service](https://earthquake.usgs.gov/fdsnws/event/1/).
- [ISC Bulletin](https://www.isc.ac.uk/iscbulletin/) and [IASPEI Seismic Format documentation](https://www.isc.ac.uk/standards/isf/).
- Wiemer, S., and Wyss, M. (2000), [Minimum magnitude of completeness in earthquake catalogs](https://doi.org/10.1785/0120000007); Aki, K. (1965), maximum-likelihood b-value estimation; Woessner, J., and Wiemer, S. (2005), magnitude completeness and catalog quality.
