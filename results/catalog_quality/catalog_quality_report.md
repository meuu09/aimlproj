# Catalog quality report

Generated from the stored catalog snapshots. This report does not modify source CSVs.
Magnitude values are grouped by their source magnitude type; no magnitude conversions are applied.
The ISC ISF depth-error values are retained with their 90% confidence metadata.
Values over 100 km and exactly 1000 km are flagged for review, not recoded or dropped.

## Catalog totals

| Catalog | Rows | IDs unique | M<3 | M≥3 | M≥4.4 | Missing depth error | >100 km depth error |
|---|---:|---:|---:|---:|---:|---:|---:|
| USGS | 566 | 566 | 0 | 566 | 299 | 286 | 0 |
| ISC v1 | 3024 | 0 | 1133 | 1676 | 268 | 1778 | 112 |
| combined v1 | 2839 | 2839 | 1133 | 1706 | 331 | 1654 | 87 |
| ISC v2 | 3024 | 0 | 1133 | 1676 | 269 | 1747 | 112 |
| combined v2 | 2839 | 2839 | 1133 | 1706 | 331 | 1559 | 99 |

## Magnitude scales

See `magnitude_type_summary.csv` for counts by catalog, source, and exact `magType`.

- USGS / USGS / mw: 4 rows, 4 at M≥3, 4 at M≥4.4.
- USGS / USGS / mb: 543 rows, 543 at M≥3, 278 at M≥4.4.
- USGS / USGS / ms: 2 rows, 2 at M≥3, 2 at M≥4.4.
- USGS / USGS / md: 1 rows, 1 at M≥3, 0 at M≥4.4.
- USGS / USGS / mwc: 3 rows, 3 at M≥3, 3 at M≥4.4.
- USGS / USGS / mww: 10 rows, 10 at M≥3, 10 at M≥4.4.
- USGS / USGS / mwb: 1 rows, 1 at M≥3, 1 at M≥4.4.
- USGS / USGS / ml: 1 rows, 1 at M≥3, 0 at M≥4.4.
- USGS / USGS / mwr: 1 rows, 1 at M≥3, 1 at M≥4.4.
- ISC v2 / ISC / type missing: 215 rows, 0 at M≥3, 0 at M≥4.4.
- ISC v2 / ISC / mb: 1002 rows, 992 at M≥3, 232 at M≥4.4.
- ISC v2 / ISC / MW: 321 rows, 221 at M≥3, 29 at M≥4.4.
- ISC v2 / ISC / MB: 2 rows, 2 at M≥3, 2 at M≥4.4.
- ISC v2 / ISC / mL: 12 rows, 7 at M≥3, 1 at M≥4.4.
- ISC v2 / ISC / MD: 555 rows, 94 at M≥3, 1 at M≥4.4.
- ISC v2 / ISC / ML: 907 rows, 350 at M≥3, 0 at M≥4.4.
- ISC v2 / ISC / mpv: 6 rows, 6 at M≥3, 0 at M≥4.4.
- ISC v2 / ISC / Mw: 4 rows, 4 at M≥3, 4 at M≥4.4.
- combined v2 / ISC / mb: 533 rows, 523 at M≥3, 30 at M≥4.4.
- combined v2 / USGS+ISC / mw: 4 rows, 4 at M≥3, 4 at M≥4.4.
- combined v2 / USGS+ISC / mb: 514 rows, 514 at M≥3, 269 at M≥4.4.
- combined v2 / USGS+ISC / ms: 2 rows, 2 at M≥3, 2 at M≥4.4.
- combined v2 / USGS / mb: 29 rows, 29 at M≥3, 9 at M≥4.4.
- combined v2 / ISC / mL: 12 rows, 7 at M≥3, 1 at M≥4.4.
- combined v2 / ISC / MD: 555 rows, 94 at M≥3, 1 at M≥4.4.
- combined v2 / ISC / ML: 906 rows, 349 at M≥3, 0 at M≥4.4.
- combined v2 / USGS+ISC / md: 1 rows, 1 at M≥3, 0 at M≥4.4.
- combined v2 / USGS+ISC / mwc: 2 rows, 2 at M≥3, 2 at M≥4.4.
- combined v2 / ISC / mpv: 6 rows, 6 at M≥3, 0 at M≥4.4.
- combined v2 / USGS / mwc: 1 rows, 1 at M≥3, 1 at M≥4.4.
- combined v2 / USGS+ISC / mww: 10 rows, 10 at M≥3, 10 at M≥4.4.
- combined v2 / USGS+ISC / mwb: 1 rows, 1 at M≥3, 1 at M≥4.4.
- combined v2 / USGS+ISC / ml: 1 rows, 1 at M≥3, 0 at M≥4.4.
- combined v2 / ISC / MW: 261 rows, 161 at M≥3, 0 at M≥4.4.
- combined v2 / USGS+ISC / mwr: 1 rows, 1 at M≥3, 1 at M≥4.4.

## Depth quality

Depth categories retain the catalog's fixed, depth-phase-constrained, solved, and unknown distinctions.
An ISC depth may be supplied by an origin other than the event's prime origin; see `depth_origin_diagnostics.csv`.
Missing error bars remain missing; they are not treated as zero uncertainty.
See `depth_origin_summary.csv` for counts by catalog of prime/depth author differences, extreme depth errors, and origin-coordinate offsets.
Depth errors equal to 1000 km remain unchanged and flagged for review; the report does not infer a missing-value sentinel.

## Event matching

The accepted one-to-one rule is ≤20 s and ≤100 km. It accepted 536 pairs; 130 have review flags, 129 exceed 20 km, and 12 exceed 50 km. Flags do not delete or alter matches.

## Filtering

See `filtering_summary.csv` for decade and 1990+ counts before and after the solved-depth case selection.
The combined catalogs contain rows below M3; analyses with an M3 floor report those rows as excluded.
