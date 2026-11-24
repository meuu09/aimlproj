# Scientific assumptions, parameter choices and their limits

Each entry gives the choice, the reason, alternatives tested, and the limitation. Numbers refer to files in `results/`.

## Data and catalog

| Choice | Reason | Alternatives / test | Limitation |
|---|---|---|---|
| Study box 31-35°N, 74-79.5°E | Most complete of four candidate boxes; continuous NW Himalayan arc segment around Kishtwar-Chamba | India-only clipping rejected: political boundaries are not geological ones. The box includes areas administered by Pakistan and China, which is stated in the README | The assignment asks for "mostly" Indian Himalaya, and most solved intermediate events lie in J&K and Himachal |
| Catalog period 1970-01-01 to 2026-11-15 (exclusive), analysis from 1990 | Download range; 1990+ for rates | Scenarios from 2000 | Completeness varies (below) |
| One whole origin per event (Task 8) | A hypocentre is a joint solution: depth trades off with origin time and epicentre | Hierarchy: USGS origin if USGS solved the depth; else first trusted ISC agency (ISC-EHB > ISC > NEIC > NDI) with a free depth and an error; else USGS. Preferring ISC-EHB over solved USGS origins was considered but not adopted, to keep established rows stable | Different agencies use different velocity models |
| Trusted depth agencies ISC-EHB, ISC, NEIC, NDI | ISC-EHB is a reprocessed catalog with depth phases (Engdahl et al. 1998); ISC and NEIC are global relocations; NDI is the national network | Other agencies' free depths are recorded (`n_free_depth_origins`) but not used | NDI errors range up to sentinel values |
| Error values ≤ 0 or ≥ 999 km set missing | NDI reports 1000 km (44 events); ISC ellipses up to 10000 km | Kept values 50-999 km (93 events) | A missing error makes an event "uncertain" |
| USGS-ISC match: 20 s, 100 km, greedy one-to-one on the ISC prime origin | Original choice; sensitivity: 10 s/50 km gives 524, 30 s/150 km gives 538 matches | 12 matches > 50 km reviewed; 2 questionable kept and flagged | Old (pre-2000) epicentres differ by tens of km between agencies |
| Duplicates: within 60 s and 50 km, only ISC-only rows flagged | Avoids merging two USGS events through a chain | `--drop-duplicates` sensitivity: results unchanged | 7 rows; rule is heuristic |
| Magnitudes kept on their own scale; `magType_norm` groups labels of the same scale | Converting mb/ML/MD to Mw needs regional relations we do not have | mb-only analyses (homogeneous scale) | The combined catalog mixes scales; its b and Mc are not physical |

## Depth

| Choice | Reason | Alternatives / test | Limitation |
|---|---|---|---|
| Intermediate = 40 ≤ depth − error and depth + error ≤ 150 km (Task 16) | Paper: ">40 km", "40-100 km", observed "40 to 110 km". The interval rule uses the error at both bounds | Lower bounds 40/50/60/70 km, upper 100/150 km (`robustness/depth_threshold_sensitivity.csv`) | **The Kishtwar band shrinks from 11 to 3 (USGS) and from 17 to 5 (combined) events at a 50 km lower bound** |
| Missing depth error means unknown, never zero (Task 18) | A missing error is not a perfect depth | Such events are "uncertain"; the Monte Carlo leaves them unperturbed and flags them | 12 USGS / 46 combined solved events lack errors |
| Reported depth and horizontal errors treated as 1σ (Monte Carlo baseline) | The USGS/ISC documentation could not be retrieved to confirm the definition | Alternatives: errors as 90 % bounds (σ = err/1.645), errors understated ×2, missing horizontal errors = 2 × the median USGS horizontal error (8.1 km) | With errors ×2, recovery of ≥90 % of the band drops to 4-15 % |
| Fixed / unflagged depths excluded in the main case | USGS `operator assigned` = analyst-fixed; unflagged events include 45 km values with zero error | Including-fixed case kept everywhere | The share of solved depths changes by decade (16/33/46/22 % at M ≥ 4.4) |
| Moho context | Sharma et al. (2018, EGU2018-11932, J&K broadband array): crustal thickness "~44 km" beneath the Lesser Himalaya, "increasing gradually … to ~54 km beneath the Higher Himalaya". Rai et al. (2006, GRL 33, L15308) give a northward-deepening Indian Moho beneath the NW Himalaya and Ladakh (values seen only in a secondary summary, not verified) | No gridded Moho model applied | Events at 40-60 km ± ~10 km cannot be assigned to the lower crust or the mantle. "Intermediate depth" here does not mean "mantle" |
| Global CMT centroid depths as an independent check | Waveform-modelled source depths | 24 solutions, all ≤ 41 km; the only matched catalog-intermediate event (49.9 ± 2.9 km) has a centroid at 25 km | Centroid ≠ hypocentre; GCMT depths of M~5 events are often weakly constrained |

## Completeness and Gutenberg-Richter

| Choice | Reason | Alternatives / test | Limitation |
|---|---|---|---|
| Authors' Mc: eqcorrscan 0.4.2 maximum curvature (0.5 bins) + 0.4 | Paper: Wiemer & Wyss (2000) maximum curvature with the Woessner & Wiemer (2005) +0.4 | Wiemer-Wyss mode (0.1 bins), GFT 90/95 % (Wiemer & Wyss 2000), MBS (Cao & Gao 2002) | GFT and MBS are mostly undefined on our data: the distributions are not clean G-R |
| M 4.4 = exploratory threshold | Authors' Mc over 1990+ = 4.0 + 0.4 | M 4.9 (USGS-complete in every decade), M 4.4 from 2000, mb 3.9 (combined mb-only, complete from the 1990s) | USGS 1990s Mc = 4.9: M 4.4 rates for the 1990s are incomplete |
| b by MLE (Aki 1965; Utsu 1965) with binning correction Δm = 0.1 | All magnitudes are on a 0.1 grid (`mc_methods.csv`) | Released least-squares fit for comparison; known bias (Marzocchi & Sandri 2003) | — |
| Report b only for N ≥ 20; "reliable" for N ≥ 50 | Aki σ = b/√N is about 22 % of b at N = 20 and 14 % at N = 50 | Shi & Bolt (1982) σ; bootstrap 95 % interval | Sampling error only; the bootstrap assumes independent events |
| b depends on the cutoff (mb b 1.2-1.3 at 4.2 to 2.1-2.2 at 4.6) | Measured | Candidates: mb saturation, incompleteness, small N, mixed scales | Cause not determined; not asserted |

## Temporal analysis

| Choice | Reason | Alternatives / test | Limitation |
|---|---|---|---|
| Rates with exact Poisson CIs (Garwood 1936) | Small counts | Annual dispersion index (variance/mean) against Poisson | Overdispersion in the all-event catalogs (aftershocks, completeness) |
| CV compared with simulated stationary Poisson CVs (same N, same window, 2000 simulations, seed 0) | A CV > 1 alone does not prove clustering | — | Changing completeness inflates the CV; intermediate events: p = 0.03-0.65 (mostly > 0.05, too few events) |
| Båth's law only within candidate sequences | Båth (1965) concerns a mainshock and its largest aftershock | The paper-style "largest minus second-largest in the cluster" is kept and labelled | No intermediate sequence: not applicable |
| Omori's law not fitted | No aftershock series at intermediate depth | — | Blocked by the data |

## Clustering and migration

| Choice | Reason | Alternatives / test | Limitation |
|---|---|---|---|
| DBSCAN (Ester et al. 1996), haversine, eps 30 km, min 4 | eps is about the location-error scale | Grid 15-45 km × 3-6 events; 3-D Euclidean DBSCAN in km (local tangent plane, < 1 % distortion over the box); k-distance plot | Clusters are density descriptions, not structures; the band's best cluster is heavily contaminated in the Monte Carlo |
| Kishtwar outline drawn by hand | Drawn around solved intermediate events (risk of circularity) | DBSCAN-independent recovery; buffered/shrunk outlines | Membership changes by ±40 % with ±0.2° buffers |
| PCA strike on km coordinates, axis 0-180° | Units comparable on both axes | Variants; GEM fault GAF_534 trace 131° as an independent reference | GEM coverage of the region is incomplete |
| Migration window 20 days / 20 km, 3-D, from any first event, minimum 3 events | Paper Sect. 3.2.3 (minimum count not stated) | 5-60 days × 10-50 km grid; all candidates kept; strict subsets dropped only after saving | Depth errors (8-15 km) are comparable to the window; the detection threshold is about M 4 |
| Linear (free and through-origin) and diffusive r = √(4πDt) fits | Shapiro et al. (1997), as cited in the paper | Permutation test of the distance-time trend | Not applied at intermediate depth (no candidates) |
| Space-time clustering (ST-DBSCAN; Birant & Kut 2007) not used | 0-1 core points within 1 year | Feasibility table | — |

## Interpretation (hypotheses, not results)

- **Metamorphic dehydration / eclogitization** (paper; Hacker et al. 2003 for dehydration embrittlement). Our data neither support nor contradict it. The catalog cannot resolve the mechanism, and our "intermediate" events may be lower-crustal.
- **Faulting near the Main Himalayan Thrust ramp and its splays.** Sharma et al. (2018) place the 2013 Kishtwar mb 5.7 "at the juncture of the MHT and the upward splay of the MBT". Shamim et al. (2024, EGU24-15623) find local-network Kishtwar seismicity mostly shallower than 25 km. GCMT mechanisms are mostly thrust (12 of 24) with shallow centroids. This explanation is consistent with the shallow part of the seismicity.
- **Depth artefact.** Bulletin depths at 40-60 km with errors of about 10 km, GCMT centroids ≤ 41 km, and the collapse of the band at a 50 km bound together leave open that some "intermediate" events are mislocated in depth. This is unresolved, and relocation with local data (as in Shamim et al. 2024) would be needed.

## References

Aki (1965) Bull. Earthq. Res. Inst. 43; Båth (1965) Tectonophysics 2; Birant & Kut (2007) Data Knowl. Eng. 60;
Cao & Gao (2002) GRL 29(9); Dziewonski et al. (1981) JGR 86; Ekström et al. (2012) PEPI 200-201;
Engdahl, van der Hilst & Buland (1998) BSSA 88; Ester et al. (1996) KDD-96; Garwood (1936) Biometrika 28;
Hacker et al. (2003) JGR 108(B1); Marzocchi & Sandri (2003) Ann. Geophys. 46; Michailos, Carpenter & Hetényi (2021)
Front. Earth Sci. 9:742700; Mignan & Woessner (2012) CORSSA; Rai et al. (2006) GRL 33 L15308;
Shamim et al. (2024) EGU24-15623; Shapiro, Huenges & Borm (1997) GJI 131; Sharma et al. (2018) EGU2018-11932;
Shi & Bolt (1982) BSSA 72; Styron & Pagani (2020) Earthq. Spectra 36; Utsu (1965) Geophys. Bull. Hokkaido Univ. 13;
Wiemer & Wyss (2000) BSSA 90; Woessner & Wiemer (2005) BSSA 95; Zoback (1992) JGR 97(B8).
