# Paper methodology vs this project

Paper: Michailos, Carpenter & Hetényi (2021), *Front. Earth Sci.* 9:742700. Quotes are from the published PDF.

| Method | What the authors did | What we do | Status | Why different | Affects interpretation? |
|---|---|---|---|---|---|
| Catalog construction | Waveforms from HIMNT, Hi-CLIMB and BPE, 2001-2003. kpick, EQTransformer and EQcorrscan template matching; NonLinLoc and HypoDD. "414 robust earthquake locations". | Bulletin catalogs: USGS ComCat and ISC (all agencies), 1970-2026, one whole origin per event | Adapted (waveform steps out of scope) | No waveforms; the brief limits us to catalogs | Yes: our depths are bulletin depths with errors of about 8-15 km, not relocations |
| Magnitudes | One local scale: "ML ranging between ML 0.8 and 4.5" | Reported magnitudes: USGS mostly mb; ISC mb/ML/MD/Mw; normalised labels; scales never merged | Adapted | No common waveform-based magnitude | Yes: b-values and Mc not comparable in absolute terms |
| Intermediate depth | ">40 km" (abstract), "40-100 km" (introduction), observed "40 to 110 km" | 40 ≤ depth ± error ≤ 150 km, plus a sensitivity analysis for lower bounds of 40-70 km and upper bounds of 100/150 km | Extended | Bulletin depths need their errors taken into account | Yes: the Kishtwar band depends strongly on the 40 km bound |
| Clusters | Cluster 1 = "box in Figure 2" (248 events); Cluster 2 near Udayapur (61) | Kishtwar outline drawn by hand; DBSCAN (2-D and 3-D) as an independent check; polygon variants | Adapted and extended | Different region; objective clustering added | Outline membership changes with ±0.1-0.2° buffers |
| Mc | "maximum curvature method (Wiemer and Wyss, 2000), applying 0.4 magnitude units upward correction" | Same (released eqcorrscan code + 0.4); also the Wiemer-Wyss mode, GFT and MBS; per catalog, decade, type and cluster | Reproduced and extended | Completeness changes over 56 years and across agencies | Yes: M 4.4 is incomplete in the 1990s (USGS 4.9) |
| b-value | "fitting the cumulative density function with a power-law above the magnitude of completeness"; 0.91 ± 0.004 (full), 0.96 (Cluster 1) | Released least-squares fit, plus MLE with Aki, Shi-Bolt and bootstrap errors; suppressed below N = 20 | Reproduced and extended | Least-squares fits to cumulative counts are biased; small N | Our b-values (≈1.2-1.4) are mixed-scale and threshold-sensitive |
| Rates and inter-event times | Fig. 5 cumulative numbers, Fig. 6 inter-event times; rates "remarkably continuous" | Same plots; rates with exact Poisson CIs; CV tested against simulated Poisson; annual dispersion | Reproduced and extended | Decades of data with changing completeness | Intermediate rates are not distinguishable from Poisson |
| Båth's law | Largest minus next-largest in the cluster: "very small (e.g., 0.2 magnitude units)" compared with Båth's 1.2 | The paper-style value is kept and labelled; a true Båth comparison is made only within candidate sequences (none at intermediate depth) | Reproduced; method corrected | Unrelated events are not a mainshock-aftershock pair | The paper-style number is not a test of Båth's law |
| Migration search | "within a 3D volume with respect to any first event (time and space window of 20 days and 20 km)"; 4 sequences; v < 1.0 km/day, D < 1.0 m²/s | Same window (minimum 3 events, not stated in the paper); all candidates kept; region-wide; window grid; free and origin linear fits, fitted and envelope D, permutation test | Reproduced and extended | Search code not released | 0 intermediate candidates in any window |
| Diffusion | Shapiro et al. 1997 cited | r = √(4πDt) (Shapiro et al. 1997); envelope and least-squares D | Reproduced | - | Not applicable without sequences |
| Figures 8-10 | Space-time plots along HypoDD cross-sections; 3-D view with focal mechanisms | Along/across-strike vs time on a PCA axis; 3-D in km with stated ×2 exaggeration; base map with GCMT mechanisms | Adapted | No relocations; PCA axis is objective | - |
| Focal mechanisms | Global CMT and earlier studies (Fig. 1); one full moment tensor with an isotropic component (Alvizuri & Hetényi 2019) | 24 Global CMT solutions with P/T/B axes and WSM regime (Zoback 1992) | Adapted | Public solutions only | All centroids ≤ 41 km: none confirms intermediate depth |
| Interpretation | Metamorphic dehydration and eclogitization; alternative dextral basement fault | Hypotheses listed against our observations; not tested | Not applicable (catalog cannot test) | - | - |

## Central vs north-western Himalaya: what can be compared

| Aspect | Paper (central) | This project (NW) | Comparable? |
|---|---|---|---|
| Period | late 2001 - mid 2003 | 1970-2026 (analysis 1990+) | No: different durations and networks |
| Magnitude | ML 0.8-4.5, Mc 2.4 | mb/ML/MD/Mw, Mc 3.9-4.9 | No: b and rates not comparable numerically |
| Locations | Relocated (HypoDD) | Bulletin, depth error ~8-15 km | Depth comparisons qualitative only |
| Intermediate events | 414 in total (248 in Cluster 1) | solved intermediate: 28 (USGS) / 54 (combined) all years; 25 / 50 since 1990 | Only qualitatively |
| Mainshock-aftershock pattern | None found | None at intermediate depth | Same qualitative result, different power |
| Migration | 4 sequences | 0 | Not comparable: detection threshold ~M4 vs ~M1 |
