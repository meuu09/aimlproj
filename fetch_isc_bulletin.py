"""Download and parse the ISC Bulletin (ISF) without altering legacy CSVs.

By default this script writes versioned outputs:
``data/isc_bulletin_events_v2.csv`` and ``data/isc_bulletin_origins_v2.csv``.
Use ``--parse-only`` to rebuild those tables from the already stored annual raw
files without making network requests.

ISF depth errors are 90% errors in kilometres. Their confidence metadata is
preserved; it is not treated as interchangeable with the USGS ``depthError``.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import time
import urllib.request
from typing import Any

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "data", "isc_raw")
OUT_EVENTS = os.path.join(HERE, "data", "isc_bulletin_events_v2.csv")
OUT_ORIGINS = os.path.join(HERE, "data", "isc_bulletin_origins_v2.csv")
BOX = dict(bot_lat=31, top_lat=35, left_lon=74, right_lon=79.5)
YEARS = range(1970, 2027)
DEPTH_PRIORITY = ["ISC-EHB", "ISC", "NEIC", "NDI"]
# Priorities are normalized to case-folded magnitude names, but the original
# type string and value are kept unchanged in the output.
MAG_PRIORITY = [
    ("mw", None),
    ("mb", "isc"),
    ("mb", "neic"),
    ("mb", None),
    ("ml", "ndi"),
    ("ml", None),
    ("ms", None),
]
URL = (
    "https://www.isc.ac.uk/cgi-bin/web-db-run?request=COMPREHENSIVE&out_format=ISF"
    "&searchshape=RECT&bot_lat={bot_lat}&top_lat={top_lat}&left_lon={left_lon}&right_lon={right_lon}"
    "&start_year={y}&start_month=1&start_day=01&start_time=00%3A00%3A00"
    "&end_year={y}&end_month=12&end_day=31&end_time=23%3A59%3A59"
    "&req_mag_type=Any&req_mag_agcy=Any&include_magnitudes=on&include_headers=on&include_comments=on"
)
ORIGIN_RE = re.compile(r"^\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}")


def fnum(value: str) -> float | None:
    value = value.strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def download(year: int) -> str:
    path = os.path.join(RAW, f"isc_isf_{year}.txt")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    for attempt in range(4):
        try:
            with urllib.request.urlopen(URL.format(y=year, **BOX), timeout=600) as response:
                content = response.read().decode("utf-8", "replace")
            if "STOP" in content or "No events were found" in content:
                with open(path, "w", encoding="utf-8", newline="") as output:
                    output.write(content)
                return path
        except Exception as exc:  # network errors are retried, never hidden
            print(f"  {year}: attempt {attempt + 1} failed ({type(exc).__name__})", flush=True)
        time.sleep(15 * (attempt + 1))
    raise RuntimeError(f"ISC download failed for {year}")


def parse(path: str):
    """Yield one event dictionary with its original ISC origin/magnitude rows."""
    event: dict[str, Any] | None = None
    last_origin: dict[str, Any] | None = None
    mode: str | None = None
    with open(path, encoding="utf-8", errors="replace") as source:
        for raw_line in source:
            line = raw_line.rstrip("\r\n")
            if line.startswith("Event "):
                if event:
                    yield event
                parts = line.split(None, 2)
                event = {
                    "id": parts[1],
                    "region": parts[2] if len(parts) > 2 else "",
                    "origins": [],
                    "mags": [],
                }
                last_origin = None
                mode = None
                continue
            if event is None:
                continue
            if line.lstrip().startswith("Date       Time"):
                mode = "origin"
                continue
            if line.startswith("Magnitude"):
                mode = "mag"
                continue
            if mode == "origin" and ORIGIN_RE.match(line):
                # ISF uses 1-based fixed columns. The depth flag is f (fixed),
                # d (fixed to a depth-phase depth), or blank (not fixed).
                depth_flag = line[76:77].strip().lower()
                origin = {
                    "time": line[0:10].replace("/", "-") + "T" + line[11:22].strip(),
                    "lat": fnum(line[36:44]),
                    "lon": fnum(line[45:54]),
                    "depth": fnum(line[71:76]),
                    "depth_flag": depth_flag,
                    "depth_fixed": depth_flag == "f",
                    "depth_err": fnum(line[78:82]),
                    "ndef": fnum(line[83:87]),
                    "author": line[118:127].strip(),
                    "origid": line[128:136].strip(),
                    "prime": False,
                    "ehb_quality": "",
                }
                event["origins"].append(origin)
                last_origin = origin
                continue
            if mode == "origin" and "(#PRIME)" in line and last_origin is not None:
                last_origin["prime"] = True
                continue
            if mode == "origin" and "DEPTH_QUALITY=" in line and last_origin is not None:
                last_origin["ehb_quality"] = line.split("DEPTH_QUALITY=")[1].strip(" )")
                continue
            if mode == "mag" and line.strip() and not line.startswith(" "):
                # Fixed columns preserve the magnitude error and station count
                # instead of dropping them while tokenizing the row.
                magnitude = fnum(line[6:10])
                if magnitude is not None:
                    event["mags"].append(
                        {
                            "type": line[0:5].strip(),
                            "mag": magnitude,
                            "error": fnum(line[11:14]),
                            "nsta": fnum(line[15:19]),
                            "author": line[20:29].strip(),
                            "origid": line[30:38].strip(),
                        }
                    )
    if event:
        yield event


def _depth_rank(origin: dict[str, Any]) -> int:
    author = str(origin.get("author", "")).upper()
    return DEPTH_PRIORITY.index(author) if author in DEPTH_PRIORITY else len(DEPTH_PRIORITY)


def _magnitude_rank(magnitude: dict[str, Any]) -> int:
    mag_type = str(magnitude.get("type", "")).casefold()
    author = str(magnitude.get("author", "")).casefold()
    for index, (wanted_type, wanted_author) in enumerate(MAG_PRIORITY):
        if mag_type == wanted_type and (wanted_author is None or author == wanted_author):
            return index
    return len(MAG_PRIORITY)


def summarise(event: dict[str, Any]) -> dict[str, Any] | None:
    origins = [origin for origin in event["origins"] if origin["lat"] is not None]
    if not origins:
        return None
    prime = next((origin for origin in origins if origin["prime"]), origins[-1])

    # Blank depth flag plus a positive 90% error is a depth not fixed by the
    # agency. A "d" flag is retained as a distinct depth-phase-constrained
    # category; "f" is an operator-fixed depth.
    candidates = [
        origin
        for origin in origins
        if origin["depth"] is not None
        and (
            (origin["depth_flag"] == "" and (origin["depth_err"] or 0) > 0)
            or origin["depth_flag"] == "d"
        )
    ]
    selected_depth = min(candidates, key=_depth_rank) if candidates else prime
    depth_flag = selected_depth["depth_flag"]
    if depth_flag == "f":
        depth_type = "operator assigned"
    elif depth_flag == "d":
        depth_type = "constrained by depth phases"
    elif selected_depth in candidates:
        depth_type = "from location"
    else:
        depth_type = "unknown"

    magnitudes = event["mags"]
    selected_magnitude = min(magnitudes, key=_magnitude_rank) if magnitudes else None
    if selected_magnitude and _magnitude_rank(selected_magnitude) == len(MAG_PRIORITY):
        # Keep the legacy fallback behavior, but preserve the exact type and
        # author so the quality report can expose its use.
        selected_magnitude = magnitudes[0]

    return {
        "isc_id": event["id"],
        "time": prime["time"],
        "latitude": prime["lat"],
        "longitude": prime["lon"],
        "prime_author": prime["author"],
        "prime_origin_id": prime["origid"],
        "depth": selected_depth["depth"],
        "depthError": selected_depth["depth_err"],
        "depth_error_confidence": "90% (ISF)",
        "depth_type": depth_type,
        "depth_author": selected_depth["author"],
        "depth_origin_id": selected_depth["origid"],
        "depth_origin_time": selected_depth["time"],
        "depth_origin_latitude": selected_depth["lat"],
        "depth_origin_longitude": selected_depth["lon"],
        "depth_flag": depth_flag,
        "depth_error_gt_100_km": bool(
            selected_depth["depth_err"] is not None and selected_depth["depth_err"] > 100
        ),
        "mag": selected_magnitude["mag"] if selected_magnitude else None,
        "magType": selected_magnitude["type"] if selected_magnitude else None,
        "mag_author": selected_magnitude["author"] if selected_magnitude else None,
        "magError": selected_magnitude["error"] if selected_magnitude else None,
        "magNsta": selected_magnitude["nsta"] if selected_magnitude else None,
        "mag_origin_id": selected_magnitude["origid"] if selected_magnitude else None,
        "n_origins": len(origins),
        "agencies": ";".join(sorted({origin["author"] for origin in origins})),
        "region": event["region"],
    }


def _write_csv(path: str, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_existing_raw_files() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    events: list[dict[str, Any]] = []
    origins: list[dict[str, Any]] = []
    for year in YEARS:
        path = os.path.join(RAW, f"isc_isf_{year}.txt")
        if not os.path.exists(path):
            continue
        for event in parse(path):
            row = summarise(event)
            if row:
                events.append(row)
            for origin in event["origins"]:
                origins.append({"isc_id": event["id"], **origin})
    events.sort(key=lambda row: (row["time"], str(row["isc_id"])))
    origins.sort(key=lambda row: (str(row["isc_id"]), row["time"], row["origid"]))
    return events, origins


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--year", type=int, help="download only this year before parsing")
    parser.add_argument(
        "--parse-only",
        action="store_true",
        help="parse existing data/isc_raw files; make no network requests",
    )
    parser.add_argument("--out-events", default=OUT_EVENTS)
    parser.add_argument("--out-origins", default=OUT_ORIGINS)
    args = parser.parse_args(argv)

    os.makedirs(RAW, exist_ok=True)
    if not args.parse_only:
        years = [args.year] if args.year is not None else list(YEARS)
        for year in years:
            download(year)
            print(f"{year}: raw bulletin available", flush=True)
            time.sleep(3)

    events, origins = parse_existing_raw_files()
    event_fields = [
        "isc_id", "time", "latitude", "longitude", "prime_author", "prime_origin_id",
        "depth", "depthError", "depth_error_confidence", "depth_type", "depth_author",
        "depth_origin_id", "depth_origin_time", "depth_origin_latitude", "depth_origin_longitude",
        "depth_flag", "depth_error_gt_100_km", "mag", "magType", "mag_author", "magError",
        "magNsta", "mag_origin_id", "n_origins", "agencies", "region",
    ]
    origin_fields = [
        "isc_id", "time", "lat", "lon", "depth", "depth_flag", "depth_fixed", "depth_err",
        "ndef", "author", "origid", "prime", "ehb_quality",
    ]
    _write_csv(args.out_events, events, event_fields)
    _write_csv(args.out_origins, origins, origin_fields)
    print(f"events={len(events)} origins={len(origins)} -> {args.out_events}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
