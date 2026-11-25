# Intermediate-depth seismicity beneath the NW (Indian) Himalaya

IIT Roorkee group assignment (topic 4). We reproduce the catalog analysis of **Michailos, Carpenter & Hetényi (2021)**, *Spatio-Temporal Evolution of Intermediate-Depth Seismicity Beneath the Himalayas*, Front. Earth Sci. 9:742700, with their code ([kemichai/Him_ides](https://github.com/kemichai/Him_ides)). We then apply it to the Himachal–Kashmir Himalaya (31–35°N, 74–79.5°E) using the open **USGS ComCat** and **ISC Bulletin** catalogs.

**Question:** do intermediate-depth earthquakes (40–150 km) here cluster in space, change in rate, occur in sequences or migrate? And how much do the catalogs limit the answer?

## Quick start

Python 3.12 is required.

```bash
uv venv --python 3.12 .venv            # or: python3.12 -m venv .venv
uv pip install --python .venv/bin/python -r requirements.txt
.venv/bin/python -m pytest -q tests     # 67 passed, 2 expected failures
.venv/bin/python run_all.py --baseline  # rebuild everything (~1 min, offline)
```

On Windows, use `py -3.12 -m venv .venv` and `.venv\Scripts\python` in place of `.venv/bin/python`.

`run_all.py` runs every step in dependency order and skips any step whose inputs failed. It logs to `results/run_log.txt` and writes `results/MANIFEST.json` (file checksums, package versions, seeds, git commit). Options:
- `--baseline`: also rerun the original code first;
- `--quick`: USGS-only analyses;
- `--download`: refresh the downloaded data.

Run a single script as a module from the repository root, for example `python -m analysis.kishtwar_cluster --catalog combined`.

## Layout

| Folder | Contents |
|---|---|
| `seis_utils.py` | Shared functions: catalog loading, depth classes, Mc, b-value, Poisson tests, matching, migration search |
| `fetch/` | Downloads: USGS catalog (new dated snapshot, never overwrites), USGS depth types, ISC Bulletin (ISF), Global CMT focal mechanisms |
| `catalog/` | USGS + ISC combined catalog, quality-control tables |
| `baseline/` | Runs the original Him_ides code and compares it with the published values |
| `analysis/` | Regional analysis, Kishtwar cluster, DBSCAN depth clusters, temporal summary, migration search, clustering checks |
| `robustness/` | Mc through time, b-value, cluster Monte Carlo and depth-threshold sensitivity |
| `tools/` | Manifest, result comparison against a git commit, before/after metrics |
| `tests/` | Unit, parser, data-consistency and baseline tests |
| `data/` | Input catalogs (raw ISC files in `data/isc_raw/`) |
| `results/`, `results_baseline/` | All generated tables and figures |
| `audit/` | Task status and before/after tables from the scientific audit |

## Data

- **USGS ComCat**, M ≥ 3, 1970–2026, with the USGS depth-type flag (solved vs analyst-fixed depth).
- **ISC Bulletin** (all agencies). Each event takes one whole origin from the first trusted agency (ISC-EHB > ISC > NEIC > NDI) that solved the depth with an error; otherwise the ISC prime origin is used.
- **Combined catalog**: USGS and ISC events matched one-to-one (≤ 20 s, ≤ 100 km). Every row is one self-consistent origin. Possible duplicates are flagged, not removed.
- **Global CMT**: 24 focal mechanisms. **GEM Global Active Faults**: subset of the box, CC BY-SA 4.0.

The box is a geological study region. It includes areas administered by Pakistan and China.

## Methods

- **Depth classes** use the error bar:
  - intermediate: depth − error ≥ 40 km and depth + error ≤ 150 km;
  - shallow: depth + error < 40 km;
  - deep: depth − error > 150 km;
  - uncertain: anything else, or no error.

  The main case uses only solved depths. A second case includes fixed depths.
- **Magnitude of completeness (Mc):** the paper's method (maximum curvature + 0.4), checked with three other methods (Wiemer–Wyss mode, goodness-of-fit test, b-value stability). M 4.4 is an exploratory cutoff, because USGS is complete only from M 4.9 in the 1990s.
- **b-value:** maximum likelihood (Aki–Utsu) with bootstrap errors; the paper's least-squares fit is kept for comparison. Not reported below 20 events.
- **Time:** rates with exact Poisson intervals. The CV of inter-event times is compared with simulated Poisson catalogs.
- **Migration:** the paper's search (20 days, 20 km, 3-D) over the whole region and a grid of windows, with linear and diffusion fits.
- **Clustering:** DBSCAN in 2-D and 3-D, with a location-error Monte Carlo (500 runs).
- Magnitudes are never converted between scales.

## Main results

1. **Baseline reproduced in part.** The full-catalog values are reproduced: 414 events, Mc 2.4, b 0.917 (published 0.91). Cluster 1's b matches at the published Mc (0.94–0.98 vs 0.96). Not reproduced:
   - the published Cluster 1 Mc of 2.85, which the released code cannot produce;
   - its count, 254 events vs 248 published;
   - the migration search, whose code was not released.
2. **Few intermediate-depth events.** Since 1990 there are 25 (USGS) and 50 (combined) with solved depths; at M ≥ 4.4 there are 8 and 15.
3. **A concentration beneath Kishtwar–Chamba** appears in the hand-drawn outline and in DBSCAN, but it is fragile:
   - it shrinks from 17 to 5 events if the lower depth bound is 50 km instead of 40;
   - it holds together in 68–79 % of Monte Carlo runs, and in only 4–15 % if the location errors are twice as large.
4. **Depth reliability is the main limitation.**
   - The crust here is about 44–54 km thick.
   - All 24 Global CMT source depths are 41 km or shallower.
   - The one catalog "intermediate" event with a CMT solution (49.9 km) has a CMT depth of 25 km.

   The catalog cannot tell lower crust from mantle.
5. **No clustering in time** is detectable for intermediate events: their inter-event times are consistent with a Poisson process. Rate changes between decades reflect completeness and the share of solved depths.
6. **No sequences or migration** at intermediate depth in any window from 5 days/10 km to 60 days/50 km. With an M ~4 detection threshold and ~10 km depth errors, this does not prove migration is absent.
7. **b-values depend on the cutoff and magnitude type**, and the distributions are not clean Gutenberg–Richter.

## Limitations

- Bulletin depth errors of about 8–15 km;
- mixed magnitude scales;
- completeness that changes over time;
- small samples;
- error definitions (1σ or otherwise) not confirmed;
- no relocation.

The paper's interpretation (metamorphic dehydration and eclogitization) cannot be tested with these catalogs.

## Licensing and attribution

- **Him_ides** has no licence, so it is not redistributed. `baseline/reproduce_baseline.py` clones it at commit 33acd38 when needed.
- **Other sources:**
  - GEM faults: CC BY-SA 4.0;
  - Global CMT (Dziewonski et al. 1981; Ekström et al. 2012);
  - ISC Bulletin (International Seismological Centre);
  - USGS ComCat (public domain).

## References

- Michailos, Carpenter & Hetényi (2021) Front. Earth Sci. 9:742700
- Aki (1965); Utsu (1965); Shi & Bolt (1982)
- Wiemer & Wyss (2000); Woessner & Wiemer (2005); Cao & Gao (2002)
- Garwood (1936); Båth (1965); Shapiro et al. (1997)
- Ester et al. (1996)
- Zoback (1992); Styron & Pagani (2020)
- Sharma et al. (2018) EGU2018-11932; Shamim et al. (2024) EGU24-15623
- ObsPy; EQcorrscan (Chamberlain et al. 2018)
