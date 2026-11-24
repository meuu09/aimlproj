"""Build a provenance-rich USGS + ISC catalog without replacing legacy data.

The default output is ``data/combined_usgs_isc_catalog_v2.csv``. Matching uses
the existing <=20 s / <=100 km, greedy one-to-one rule. Every accepted pair is
written to an audit CSV with residual and candidate-count diagnostics. Large
residuals are flagged for review and are not automatically rejected.
"""

from __future__ import annotations

import argparse
import os
from typing import Any

import numpy as np
import pandas as pd

from catalog_qc import candidate_review_reason, find_matches

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
USGS_DEFAULT = os.path.join(DATA, "usgs_himachal_kashmir_M3plus_1970-2026.csv")
DEPTH_DEFAULT = os.path.join(DATA, "usgs_depth_type.csv")
DEPTH_V2_DEFAULT = os.path.join(DATA, "usgs_depth_type_v2.csv")
ISC_V2_DEFAULT = os.path.join(DATA, "isc_bulletin_events_v2.csv")
ISC_LEGACY = os.path.join(DATA, "isc_bulletin_events.csv")
OUTPUT_DEFAULT = os.path.join(DATA, "combined_usgs_isc_catalog_v2.csv")
MATCH_DIAGNOSTICS_DEFAULT = os.path.join(HERE, "results", "catalog_quality", "event_match_diagnostics.csv")
MATCH_S, MATCH_KM = 20.0, 100.0
REVIEW_DISTANCE_KM = 20.0  # review flag only; not a rejection criterion
TRUSTED_DEPTH = ("ISC-EHB", "ISC", "NEIC", "NDI")
SOLVED = ("from location", "constrained by depth phases")

OUTPUT_COLUMNS = [
    "id", "time", "latitude", "longitude", "depth", "depthError", "depth_error_confidence",
    "depth_type", "depth_source", "depth_origin_id", "depth_origin_time", "depth_origin_latitude",
    "depth_origin_longitude", "depth_flag", "horizontalError", "mag", "magType", "magError",
    "mag_source", "magNsta", "source", "isc_id", "place", "isc_time", "isc_latitude",
    "isc_longitude", "isc_prime_author", "isc_prime_origin_id", "isc_depth", "isc_depthError",
    "isc_depth_type", "isc_depth_author", "isc_depth_origin_id", "isc_depth_origin_time",
    "isc_depth_origin_latitude", "isc_depth_origin_longitude", "isc_mag", "isc_magType",
    "isc_magError", "isc_mag_author", "isc_mag_origin_id",
]


def _is_missing(value: Any) -> bool:
    return value is None or pd.isna(value) or str(value).strip() == ""


def _float_or_none(value: Any) -> float | None:
    return None if _is_missing(value) else float(value)


def _isc_provenance(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "isc_time": row.get("time"),
        "isc_latitude": row.get("latitude"),
        "isc_longitude": row.get("longitude"),
        "isc_prime_author": row.get("prime_author"),
        "isc_prime_origin_id": row.get("prime_origin_id"),
        "isc_depth": row.get("depth"),
        "isc_depthError": row.get("depthError"),
        "isc_depth_type": row.get("depth_type"),
        "isc_depth_author": row.get("depth_author"),
        "isc_depth_origin_id": row.get("depth_origin_id"),
        "isc_depth_origin_time": row.get("depth_origin_time"),
        "isc_depth_origin_latitude": row.get("depth_origin_latitude"),
        "isc_depth_origin_longitude": row.get("depth_origin_longitude"),
        "isc_mag": row.get("mag"),
        "isc_magType": row.get("magType"),
        "isc_magError": row.get("magError"),
        "isc_mag_author": row.get("mag_author"),
        "isc_mag_origin_id": row.get("mag_origin_id"),
    }


def _base_usgs(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "id": row.get("id"),
        "time": row.get("time"),
        "latitude": row.get("latitude"),
        "longitude": row.get("longitude"),
        "depth": row.get("depth"),
        "depthError": row.get("depthError"),
        "depth_error_confidence": "USGS ComCat field; confidence level not recorded in this export",
        "depth_type": row.get("depth_type", "unknown"),
        "depth_source": "USGS",
        "depth_origin_id": None,
        "depth_origin_time": row.get("time"),
        "depth_origin_latitude": row.get("latitude"),
        "depth_origin_longitude": row.get("longitude"),
        "depth_flag": row.get("depth_type", "unknown"),
        "horizontalError": row.get("horizontalError"),
        "mag": row.get("mag"),
        "magType": row.get("magType"),
        "magError": row.get("magError"),
        "mag_source": "USGS",
        "magNsta": row.get("magNst"),
        "source": "USGS",
        "isc_id": None,
        "place": row.get("place"),
    }


def build_catalog(
    usgs_path: str,
    depth_types_path: str,
    isc_path: str,
    output_path: str,
    match_diagnostics_path: str,
) -> tuple[pd.DataFrame, list[dict[str, Any]], int]:
    usgs = pd.read_csv(usgs_path)
    depth_types = pd.read_csv(depth_types_path)[["id", "depth_type"]]
    usgs = usgs.merge(depth_types, on="id", how="left", suffixes=("", "_flag"))
    if "depth_type_flag" in usgs:
        usgs["depth_type"] = usgs["depth_type_flag"].fillna(usgs.get("depth_type", "unknown"))
        usgs = usgs.drop(columns="depth_type_flag")
    usgs["depth_type"] = usgs["depth_type"].fillna("unknown")

    isc = pd.read_csv(isc_path, dtype={"isc_id": str})
    isc = isc[isc.latitude.between(31, 35) & isc.longitude.between(74, 79.5)].reset_index(drop=True)
    isc["depth_quality_note"] = ""
    untrusted = (isc.depth_type.isin(SOLVED)) & ~isc.depth_author.isin(TRUSTED_DEPTH)
    isc.loc[untrusted, "depth_quality_note"] = "untrusted_depth_agency; original values preserved"
    isc.loc[untrusted, "depth_type"] = "unknown"

    us_records = usgs.to_dict("records")
    isc_records = isc.to_dict("records")
    matches, candidate_pair_count = find_matches(
        us_records, isc_records, max_seconds=MATCH_S, max_distance_km=MATCH_KM
    )
    us_to_match = {match["usgs_index"]: match for match in matches}
    isc_matched = {match["isc_index"] for match in matches}

    rows: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for usgs_index, us_row in enumerate(us_records):
        out = _base_usgs(us_row)
        match = us_to_match.get(usgs_index)
        if match is not None:
            isc_row = isc_records[match["isc_index"]]
            out.update(
                isc_id=isc_row.get("isc_id"),
                source="USGS+ISC",
                **_isc_provenance(isc_row),
            )
            usgs_error = _float_or_none(us_row.get("depthError"))
            usgs_solved = us_row.get("depth_type") in SOLVED and usgs_error is not None and usgs_error > 0
            if not usgs_solved and isc_row.get("depth_type") in SOLVED:
                out.update(
                    depth=isc_row.get("depth"),
                    depthError=isc_row.get("depthError"),
                    depth_error_confidence=isc_row.get("depth_error_confidence", "90% (ISF)"),
                    depth_type=isc_row.get("depth_type"),
                    depth_source=f"ISC:{isc_row.get('depth_author', '')}",
                    depth_origin_id=isc_row.get("depth_origin_id"),
                    depth_origin_time=isc_row.get("depth_origin_time"),
                    depth_origin_latitude=isc_row.get("depth_origin_latitude"),
                    depth_origin_longitude=isc_row.get("depth_origin_longitude"),
                    depth_flag=isc_row.get("depth_flag"),
                )
            reason = candidate_review_reason(match, review_distance_km=REVIEW_DISTANCE_KM)
            diagnostics.append(
                {
                    "usgs_id": us_row.get("id"),
                    "isc_id": isc_row.get("isc_id"),
                    "time_difference_seconds": match["time_difference_seconds"],
                    "horizontal_distance_km": match["horizontal_distance_km"],
                    "candidate_score": match["candidate_score"],
                    "usgs_candidate_count": match["usgs_candidate_count"],
                    "isc_candidate_count": match["isc_candidate_count"],
                    "review_required": bool(reason),
                    "review_reason": reason,
                    "usgs_mag": us_row.get("mag"),
                    "usgs_magType": us_row.get("magType"),
                    "isc_mag": isc_row.get("mag"),
                    "isc_magType": isc_row.get("magType"),
                    "isc_mag_author": isc_row.get("mag_author"),
                }
            )
        rows.append(out)

    for isc_index, isc_row in enumerate(isc_records):
        if isc_index in isc_matched or _is_missing(isc_row.get("mag")):
            continue
        out = {
            "id": f"isc{isc_row.get('isc_id')}",
            "time": isc_row.get("time"),
            "latitude": isc_row.get("latitude"),
            "longitude": isc_row.get("longitude"),
            "depth": isc_row.get("depth"),
            "depthError": isc_row.get("depthError"),
            "depth_error_confidence": isc_row.get("depth_error_confidence", "90% (ISF)"),
            "depth_type": isc_row.get("depth_type"),
            "depth_source": f"ISC:{isc_row.get('depth_author', '')}",
            "depth_origin_id": isc_row.get("depth_origin_id"),
            "depth_origin_time": isc_row.get("depth_origin_time"),
            "depth_origin_latitude": isc_row.get("depth_origin_latitude"),
            "depth_origin_longitude": isc_row.get("depth_origin_longitude"),
            "depth_flag": isc_row.get("depth_flag"),
            "horizontalError": None,
            "mag": isc_row.get("mag"),
            "magType": isc_row.get("magType"),
            "magError": isc_row.get("magError"),
            "mag_source": f"ISC:{isc_row.get('mag_author', '')}",
            "magNsta": isc_row.get("magNsta"),
            "source": "ISC",
            "isc_id": isc_row.get("isc_id"),
            "place": isc_row.get("region"),
            **_isc_provenance(isc_row),
        }
        rows.append(out)

    out = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)
    out["_sort_time"] = pd.to_datetime(out.time, utc=True, format="ISO8601")
    out = out.sort_values("_sort_time").drop(columns="_sort_time").reset_index(drop=True)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    out.to_csv(output_path, index=False)
    os.makedirs(os.path.dirname(match_diagnostics_path), exist_ok=True)
    pd.DataFrame(diagnostics).to_csv(match_diagnostics_path, index=False)
    return out, diagnostics, candidate_pair_count


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--usgs", default=USGS_DEFAULT)
    default_depth = DEPTH_V2_DEFAULT if os.path.exists(DEPTH_V2_DEFAULT) else DEPTH_DEFAULT
    parser.add_argument("--depth-types", default=default_depth)
    parser.add_argument("--isc-events", default=None)
    parser.add_argument("--output", default=OUTPUT_DEFAULT)
    parser.add_argument("--match-diagnostics", default=MATCH_DIAGNOSTICS_DEFAULT)
    args = parser.parse_args(argv)

    isc_path = args.isc_events or (ISC_V2_DEFAULT if os.path.exists(ISC_V2_DEFAULT) else ISC_LEGACY)
    output, diagnostics, candidate_count = build_catalog(
        args.usgs, args.depth_types, isc_path, args.output, args.match_diagnostics
    )
    source_counts = output.source.value_counts().to_dict()
    solved = output.depth_type.isin(SOLVED)
    print(
        f"USGS {sum(source_counts.get(k, 0) for k in ('USGS', 'USGS+ISC'))} | "
        f"ISC input {len(pd.read_csv(isc_path))} | candidate pairs {candidate_count} | "
        f"matched {len(diagnostics)} | combined {len(output)}"
    )
    print("source counts:", source_counts)
    print("depth types:", output.depth_type.value_counts(dropna=False).to_dict())
    print("depth sources (classified solved):", output.loc[solved, "depth_source"].value_counts().to_dict())
    print("M<3 rows:", int((pd.to_numeric(output.mag, errors="coerce") < 3).sum()))
    print("match diagnostics:", args.match_diagnostics)
    print("combined catalog:", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
