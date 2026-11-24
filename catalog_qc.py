"""Small, dependency-light helpers for catalog quality checks.

The matching routine deliberately preserves the existing project's candidate
limits and greedy one-to-one scoring rule. It adds diagnostics; it does not
decide that a scientifically questionable match is false.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from datetime import datetime, timedelta, timezone
from math import asin, cos, pi, sin, sqrt
from typing import Any, Iterable


EARTH_RADIUS_KM = 6371.0


def parse_utc(value: str) -> datetime:
    """Parse an ISO timestamp as UTC and truncate to the millisecond precision
    used by the original pandas matching code.
    """
    value = value.strip()
    if value.endswith("Z"):
        value = value[:-1] + "+00:00"
    result = datetime.fromisoformat(value)
    if result.tzinfo is None:
        result = result.replace(tzinfo=timezone.utc)
    else:
        result = result.astimezone(timezone.utc)
    result = result.replace(microsecond=(result.microsecond // 1000) * 1000)
    return result


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Great-circle distance in km."""
    p1, p2 = lat1 * pi / 180.0, lat2 * pi / 180.0
    dp = (lat2 - lat1) * pi / 360.0
    dl = (lon2 - lon1) * pi / 360.0
    a = sin(dp) ** 2 + cos(p1) * cos(p2) * sin(dl) ** 2
    return EARTH_RADIUS_KM * 2.0 * asin(sqrt(min(1.0, a)))


def find_matches(
    usgs_rows: list[dict[str, Any]],
    isc_rows: list[dict[str, Any]],
    *,
    max_seconds: float = 20.0,
    max_distance_km: float = 100.0,
) -> tuple[list[dict[str, Any]], int]:
    """Return the existing greedy one-to-one event matches and candidate count.

    Candidates are scored as ``abs(dt)/max_seconds + distance/max_distance``
    and sorted by (score, USGS input index, ISC input index), matching the
    original implementation's tie-breaking order.
    """
    isc_times = [parse_utc(row["time"]) for row in isc_rows]
    time_index = sorted((when, index) for index, when in enumerate(isc_times))
    sorted_times = [item[0] for item in time_index]
    window = timedelta(seconds=max_seconds)
    candidates: list[tuple[float, int, int, float, float]] = []
    usgs_degree: dict[int, int] = {}
    isc_degree: dict[int, int] = {}

    for usgs_index, usgs in enumerate(usgs_rows):
        usgs_time = parse_utc(usgs["time"])
        left = bisect_left(sorted_times, usgs_time - window)
        right = bisect_right(sorted_times, usgs_time + window)
        for isc_time, isc_index in time_index[left:right]:
            delta_seconds = abs((isc_time - usgs_time).total_seconds())
            isc = isc_rows[isc_index]
            distance = haversine_km(
                float(usgs["latitude"]),
                float(usgs["longitude"]),
                float(isc["latitude"]),
                float(isc["longitude"]),
            )
            if distance > max_distance_km:
                continue
            score = delta_seconds / max_seconds + distance / max_distance_km
            candidates.append((score, usgs_index, isc_index, delta_seconds, distance))
            usgs_degree[usgs_index] = usgs_degree.get(usgs_index, 0) + 1
            isc_degree[isc_index] = isc_degree.get(isc_index, 0) + 1

    matched_usgs: set[int] = set()
    matched_isc: set[int] = set()
    matches: list[dict[str, Any]] = []
    for score, usgs_index, isc_index, delta_seconds, distance in sorted(candidates):
        if usgs_index in matched_usgs or isc_index in matched_isc:
            continue
        matched_usgs.add(usgs_index)
        matched_isc.add(isc_index)
        matches.append(
            {
                "usgs_index": usgs_index,
                "isc_index": isc_index,
                "candidate_score": score,
                "time_difference_seconds": delta_seconds,
                "horizontal_distance_km": distance,
                "usgs_candidate_count": usgs_degree[usgs_index],
                "isc_candidate_count": isc_degree[isc_index],
            }
        )
    return matches, len(candidates)


def depth_class(
    depth: float | None,
    depth_error: float | None,
    *,
    minimum_km: float = 40.0,
    maximum_km: float = 150.0,
) -> str:
    """Apply the project's current central-depth/low-side-error convention."""
    if depth is None:
        return "uncertain"
    if depth_error is None:
        return "uncertain"
    if depth - depth_error >= minimum_km and depth <= maximum_km:
        return "intermediate"
    if depth + depth_error < minimum_km:
        return "shallow"
    return "uncertain"


def candidate_review_reason(match: dict[str, Any], *, review_distance_km: float = 20.0) -> str:
    """Return advisory diagnostics; matches remain accepted by the caller."""
    reasons: list[str] = []
    if float(match["horizontal_distance_km"]) > review_distance_km:
        reasons.append("horizontal_residual_over_review_distance")
    if int(match["usgs_candidate_count"]) > 1:
        reasons.append("multiple_isc_candidates_for_usgs")
    if int(match["isc_candidate_count"]) > 1:
        reasons.append("multiple_usgs_candidates_for_isc")
    return ";".join(reasons)


def nonempty(value: Any) -> bool:
    return value is not None and str(value).strip() != ""

