"""Write non-destructive quality and filtering summaries for stored catalogs."""

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter, defaultdict
from datetime import datetime
from statistics import median
from typing import Any

from catalog_qc import haversine_km, nonempty, parse_utc

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT_DEFAULT = os.path.join(HERE, "results", "catalog_quality")
SOLVED = {"from location", "constrained by depth phases"}


def read_rows(path: str) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8-sig") as source:
        return list(csv.DictReader(source))


def _number(value: Any) -> float | None:
    if not nonempty(value):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _year(row: dict[str, Any]) -> int | None:
    if not nonempty(row.get("time")):
        return None
    try:
        return parse_utc(str(row["time"])).year
    except (ValueError, TypeError):
        return None


def _summary_row(name: str, rows: list[dict[str, Any]]) -> dict[str, Any]:
    magnitudes = [_number(row.get("mag")) for row in rows]
    magnitudes = [value for value in magnitudes if value is not None]
    years = [_year(row) for row in rows]
    years = [value for value in years if value is not None]
    depth_errors = [_number(row.get("depthError")) for row in rows]
    mags = [_number(row.get("magError")) for row in rows]
    horizontal = [_number(row.get("horizontalError")) for row in rows]
    ids = [str(row.get("id", "")) for row in rows if nonempty(row.get("id"))]
    return {
        "catalog": name,
        "rows": len(rows),
        "unique_ids": len(set(ids)),
        "duplicate_id_rows": len(ids) - len(set(ids)),
        "first_year": min(years) if years else "",
        "last_year": max(years) if years else "",
        "magnitude_missing": len(rows) - len(magnitudes),
        "magnitude_min": min(magnitudes) if magnitudes else "",
        "magnitude_max": max(magnitudes) if magnitudes else "",
        "magnitude_lt_3": sum(value < 3.0 for value in magnitudes),
        "magnitude_ge_3": sum(value >= 3.0 for value in magnitudes),
        "magnitude_ge_4_4": sum(value >= 4.4 for value in magnitudes),
        "depth_error_missing": sum(value is None for value in depth_errors),
        "depth_error_zero": sum(value == 0 for value in depth_errors if value is not None),
        "depth_error_gt_100_km": sum(value > 100 for value in depth_errors if value is not None),
        "depth_error_eq_1000_km": sum(value == 1000 for value in depth_errors if value is not None),
        "horizontal_error_missing": sum(value is None for value in horizontal),
        "magnitude_error_missing": sum(value is None for value in mags),
    }


def _group_rows(rows: list[dict[str, Any]], keys: tuple[str, ...]):
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        values = []
        for key in keys:
            value = row.get(key)
            if not nonempty(value) and key == "depth_source":
                value = row.get("depth_author")
            if not nonempty(value) and key == "mag_source":
                value = row.get("mag_author", row.get("magSource"))
            values.append(str(value or ""))
        groups[tuple(values)].append(row)
    return groups


def build_reports(data_dir: str = DATA, output_dir: str = OUT_DEFAULT) -> dict[str, str]:
    os.makedirs(output_dir, exist_ok=True)
    paths = {
        "USGS": os.path.join(data_dir, "usgs_himachal_kashmir_M3plus_1970-2026.csv"),
        "ISC v1": os.path.join(data_dir, "isc_bulletin_events.csv"),
        "combined v1": os.path.join(data_dir, "combined_usgs_isc_catalog.csv"),
        "ISC v2": os.path.join(data_dir, "isc_bulletin_events_v2.csv"),
        "combined v2": os.path.join(data_dir, "combined_usgs_isc_catalog_v2.csv"),
    }
    catalogs: dict[str, list[dict[str, Any]]] = {}
    for name, path in paths.items():
        if os.path.exists(path):
            catalogs[name] = read_rows(path)

    # Attach USGS depth quality flags by ID without changing the source CSV.
    depth_path = os.path.join(data_dir, "usgs_depth_type.csv")
    if "USGS" in catalogs and os.path.exists(depth_path):
        flags = {row["id"]: row.get("depth_type", "unknown") for row in read_rows(depth_path)}
        for row in catalogs["USGS"]:
            row["depth_type"] = flags.get(row.get("id", ""), "unknown")
            row["source"] = "USGS"
            row["mag_source"] = row.get("magSource", "USGS") or "USGS"
    for name in ("ISC v1", "ISC v2"):
        if name in catalogs:
            for row in catalogs[name]:
                row["source"] = "ISC"

    quality = [_summary_row(name, rows) for name, rows in catalogs.items()]

    mag_rows: list[dict[str, Any]] = []
    depth_rows: list[dict[str, Any]] = []
    filter_rows: list[dict[str, Any]] = []
    for name, rows in catalogs.items():
        for (source, mag_type), subset in _group_rows(rows, ("source", "magType")).items():
            values = [_number(row.get("mag")) for row in subset]
            values = [value for value in values if value is not None]
            mag_rows.append(
                {
                    "catalog": name,
                    "source": source,
                    "magType": mag_type,
                    "rows": len(subset),
                    "magnitude_missing": len(subset) - len(values),
                    "magnitude_lt_3": sum(value < 3 for value in values),
                    "magnitude_ge_3": sum(value >= 3 for value in values),
                    "magnitude_ge_4_4": sum(value >= 4.4 for value in values),
                    "magnitude_min": min(values) if values else "",
                    "magnitude_max": max(values) if values else "",
                }
            )

        for (source, depth_type, depth_source), subset in _group_rows(
            rows, ("source", "depth_type", "depth_source")
        ).items():
            errors = [_number(row.get("depthError")) for row in subset]
            depth_rows.append(
                {
                    "catalog": name,
                    "source": source,
                    "depth_type": depth_type,
                    "depth_source": depth_source,
                    "rows": len(subset),
                    "depth_error_missing": sum(value is None for value in errors),
                    "depth_error_zero": sum(value == 0 for value in errors if value is not None),
                    "depth_error_gt_100_km": sum(value > 100 for value in errors if value is not None),
                    "depth_error_eq_1000_km": sum(value == 1000 for value in errors if value is not None),
                    "depth_error_confidence": ";".join(
                        sorted({str(row.get("depth_error_confidence", "")) for row in subset if nonempty(row.get("depth_error_confidence"))})
                    ),
                }
            )

        years_to_rows: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            year = _year(row)
            if year is not None:
                decade = f"{year // 10 * 10}s"
                years_to_rows[decade].append(row)
                if year >= 1990:
                    years_to_rows["1990+"].append(row)
        for period, subset in sorted(years_to_rows.items()):
            solved = [row for row in subset if row.get("depth_type", "unknown") in SOLVED]
            for depth_case, selected in (("all depth categories", subset), ("solved type only", solved)):
                values = [_number(row.get("mag")) for row in selected]
                values = [value for value in values if value is not None]
                filter_rows.append(
                    {
                        "catalog": name,
                        "period": period,
                        "depth_case": depth_case,
                        "rows": len(selected),
                        "magnitude_lt_3": sum(value < 3 for value in values),
                        "magnitude_ge_3": sum(value >= 3 for value in values),
                        "magnitude_ge_4_4": sum(value >= 4.4 for value in values),
                        "depth_error_missing": sum(_number(row.get("depthError")) is None for row in selected),
                    }
                )

    depth_origin_rows: list[dict[str, Any]] = []
    for name in ("ISC v2", "combined v2"):
        for row in catalogs.get(name, []):
            prime_lat, prime_lon = _number(row.get("latitude")), _number(row.get("longitude"))
            depth_lat = _number(row.get("depth_origin_latitude"))
            depth_lon = _number(row.get("depth_origin_longitude"))
            if prime_lat is None or prime_lon is None or depth_lat is None or depth_lon is None:
                distance = ""
            else:
                distance = haversine_km(prime_lat, prime_lon, depth_lat, depth_lon)
            prime_time, depth_time = row.get("time"), row.get("depth_origin_time")
            try:
                dt_seconds = abs((parse_utc(str(prime_time)) - parse_utc(str(depth_time))).total_seconds())
            except (TypeError, ValueError):
                dt_seconds = ""
            depth_origin_rows.append(
                {
                    "catalog": name,
                    "event_id": row.get("id", row.get("isc_id", "")),
                    "isc_id": row.get("isc_id", ""),
                    "depth_source": row.get("depth_source", row.get("depth_author", "")),
                    "prime_author": row.get("isc_prime_author", row.get("prime_author", "")),
                    "depth_author": row.get("isc_depth_author", row.get("depth_author", "")),
                    "prime_origin_id": row.get("isc_prime_origin_id", row.get("prime_origin_id", "")),
                    "depth_origin_id": row.get("depth_origin_id", row.get("isc_depth_origin_id", "")),
                    "prime_to_depth_origin_km": distance,
                    "prime_to_depth_origin_seconds": dt_seconds,
                    "depth_error_confidence": row.get("depth_error_confidence", ""),
                    "depthError": row.get("depthError", ""),
                }
            )

    match_path = os.path.join(output_dir, "event_match_diagnostics.csv")
    match_summary: list[dict[str, Any]] = []
    if os.path.exists(match_path):
        match_rows = read_rows(match_path)
        distances = [_number(row.get("horizontal_distance_km")) for row in match_rows]
        distances = [value for value in distances if value is not None]
        time_errors = [_number(row.get("time_difference_seconds")) for row in match_rows]
        time_errors = [value for value in time_errors if value is not None]
        match_summary.append(
            {
                "rule_time_seconds_max": 20,
                "rule_horizontal_km_max": 100,
                "accepted_matches": len(match_rows),
                "review_flagged_matches": sum(row.get("review_required", "").casefold() == "true" for row in match_rows),
                "matches_over_20_km": sum(value > 20 for value in distances),
                "matches_over_50_km": sum(value > 50 for value in distances),
                "time_residual_max_seconds": max(time_errors) if time_errors else "",
                "horizontal_residual_median_km": median(distances) if distances else "",
                "horizontal_residual_max_km": max(distances) if distances else "",
                "multiple_candidate_matches": sum(
                    (_number(row.get("usgs_candidate_count")) or 0) > 1
                    or (_number(row.get("isc_candidate_count")) or 0) > 1
                    for row in match_rows
                ),
            }
        )

    _write_csv(os.path.join(output_dir, "catalog_quality_summary.csv"), quality)
    _write_csv(os.path.join(output_dir, "magnitude_type_summary.csv"), mag_rows)
    _write_csv(os.path.join(output_dir, "depth_quality_summary.csv"), depth_rows)
    _write_csv(os.path.join(output_dir, "filtering_summary.csv"), filter_rows)
    _write_csv(os.path.join(output_dir, "depth_origin_diagnostics.csv"), depth_origin_rows)
    _write_csv(os.path.join(output_dir, "event_match_summary.csv"), match_summary)
    origin_summary = []
    for catalog_name in ("ISC v2", "combined v2"):
        subset = [row for row in depth_origin_rows if row["catalog"] == catalog_name]
        origin_summary.append(
            {
                "catalog": catalog_name,
                "rows": len(subset),
                "different_origin_author": sum(row.get("prime_author") != row.get("depth_author") for row in subset),
                "depthError_eq_1000_km": sum(_number(row.get("depthError")) == 1000 for row in subset),
                "prime_depth_origin_over_20_km": sum(
                    (_number(row.get("prime_to_depth_origin_km")) or 0) > 20 for row in subset
                ),
            }
        )
    _write_csv(os.path.join(output_dir, "depth_origin_summary.csv"), origin_summary)
    _write_markdown(output_dir, quality, mag_rows, depth_rows, match_summary, origin_summary)
    return {name: os.path.join(output_dir, filename) for name, filename in {
        "catalog_quality": "catalog_quality_summary.csv",
        "magnitude_types": "magnitude_type_summary.csv",
        "depth_quality": "depth_quality_summary.csv",
        "filtering": "filtering_summary.csv",
        "depth_origins": "depth_origin_diagnostics.csv",
        "depth_origin_summary": "depth_origin_summary.csv",
        "match_summary": "event_match_summary.csv",
        "report": "catalog_quality_report.md",
    }.items()}


def _write_csv(path: str, rows: list[dict[str, Any]]) -> None:
    if not rows:
        with open(path, "w", newline="", encoding="utf-8") as output:
            output.write("\n")
        return
    columns = list(rows[0])
    with open(path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_markdown(
    output_dir: str,
    quality: list[dict[str, Any]],
    mag_rows: list[dict[str, Any]],
    depth_rows: list[dict[str, Any]],
    match_summary: list[dict[str, Any]],
    origin_summary: list[dict[str, Any]],
) -> None:
    lines = [
        "# Catalog quality report",
        "",
        "Generated from the stored catalog snapshots. This report does not modify source CSVs.",
        "Magnitude values are grouped by their source magnitude type; no magnitude conversions are applied.",
        "The ISC ISF depth-error values are retained with their 90% confidence metadata.",
        "Values over 100 km and exactly 1000 km are flagged for review, not recoded or dropped.",
        "",
        "## Catalog totals",
        "",
        "| Catalog | Rows | IDs unique | M<3 | M≥3 | M≥4.4 | Missing depth error | >100 km depth error |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in quality:
        lines.append(
            f"| {row['catalog']} | {row['rows']} | {row['unique_ids']} | {row['magnitude_lt_3']} | "
            f"{row['magnitude_ge_3']} | {row['magnitude_ge_4_4']} | {row['depth_error_missing']} | "
            f"{row['depth_error_gt_100_km']} |"
        )
    lines.extend(["", "## Magnitude scales", "", "See `magnitude_type_summary.csv` for counts by catalog, source, and exact `magType`.", ""])
    for row in mag_rows:
        if row["catalog"] in {"USGS", "ISC v2", "combined v2"}:
            lines.append(
                f"- {row['catalog']} / {row['source'] or 'source untagged'} / {row['magType'] or 'type missing'}: "
                f"{row['rows']} rows, {row['magnitude_ge_3']} at M≥3, {row['magnitude_ge_4_4']} at M≥4.4."
            )
    lines.extend(
        [
            "",
            "## Depth quality",
            "",
            "Depth categories retain the catalog's fixed, depth-phase-constrained, solved, and unknown distinctions.",
            "An ISC depth may be supplied by an origin other than the event's prime origin; see `depth_origin_diagnostics.csv`.",
            "Missing error bars remain missing; they are not treated as zero uncertainty.",
            "See `depth_origin_summary.csv` for counts by catalog of prime/depth author differences, extreme depth errors, and origin-coordinate offsets.",
            "Depth errors equal to 1000 km remain unchanged and flagged for review; the report does not infer a missing-value sentinel.",
            "",
            "## Event matching",
            "",
        ]
    )
    if match_summary:
        match = match_summary[0]
        lines.append(
            f"The accepted one-to-one rule is ≤20 s and ≤100 km. It accepted {match['accepted_matches']} pairs; "
            f"{match['review_flagged_matches']} have review flags, {match['matches_over_20_km']} exceed 20 km, "
            f"and {match['matches_over_50_km']} exceed 50 km. Flags do not delete or alter matches."
        )
    else:
        lines.append("Run `build_combined_catalog.py` to create match residual diagnostics.")
    lines.extend(
        [
            "",
            "## Filtering",
            "",
            "See `filtering_summary.csv` for decade and 1990+ counts before and after the solved-depth case selection.",
            "The combined catalogs contain rows below M3; analyses with an M3 floor report those rows as excluded.",
            "",
        ]
    )
    with open(os.path.join(output_dir, "catalog_quality_report.md"), "w", encoding="utf-8") as output:
        output.write("\n".join(lines))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", default=DATA)
    parser.add_argument("--output-dir", default=OUT_DEFAULT)
    args = parser.parse_args(argv)
    outputs = build_reports(args.data_dir, args.output_dir)
    for name, path in outputs.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
