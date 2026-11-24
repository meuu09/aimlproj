"""Fetch the USGS ComCat origin depth-type flag for each stored event.

Defaults use versioned inputs/outputs so the original snapshot and flag table
are not replaced. Failed requests are retried on the next run.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.request
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SNAPSHOT_V2 = os.path.join(DATA, "usgs_himachal_kashmir_M3plus_1970-2026_snapshot_v2.csv")
SOURCE_LEGACY = os.path.join(DATA, "usgs_himachal_kashmir_M3plus_1970-2026.csv")
OUTPUT_DEFAULT = os.path.join(DATA, "usgs_depth_type_v2.csv")
URL = "https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={}&format=geojson"
FIELDS = ["id", "depth_type", "origin_source", "depth_km", "origin_product_id", "retrieved_utc"]


def fetch_flags(source_path: str, output_path: str, *, pause_seconds: float = 0.2) -> int:
    done: dict[str, dict[str, str]] = {}
    if os.path.exists(output_path):
        with open(output_path, newline="", encoding="utf-8-sig") as existing:
            for row in csv.DictReader(existing):
                if row.get("depth_type") != "fetch_failed":
                    done[row["id"]] = row
    with open(source_path, newline="", encoding="utf-8-sig") as source:
        rows = list(csv.DictReader(source))

    write_header = not os.path.exists(output_path) or os.path.getsize(output_path) == 0
    with open(output_path, "a", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=FIELDS, extrasaction="ignore")
        if write_header:
            writer.writeheader()
        for index, row in enumerate(rows, start=1):
            event_id = row["id"]
            if event_id in done:
                continue
            result = None
            for attempt in range(3):
                try:
                    request = urllib.request.Request(
                        URL.format(event_id), headers={"User-Agent": "himalayan-seismicity-research/1.0"}
                    )
                    with urllib.request.urlopen(request, timeout=30) as response:
                        payload = json.load(response)
                    products = payload.get("properties", {}).get("products", {}).get("origin", [])
                    if not products:
                        raise ValueError("ComCat response has no origin product")
                    origin = products[0]
                    properties = origin.get("properties", {})
                    depth_type = str(properties.get("depth-type") or "unknown").lower().replace("_", " ")
                    result = {
                        "id": event_id,
                        "depth_type": depth_type,
                        "origin_source": properties.get("origin-source", ""),
                        "depth_km": properties.get("depth", ""),
                        "origin_product_id": origin.get("id", ""),
                        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                    }
                    break
                except Exception as exc:
                    print(f"{event_id}: attempt {attempt + 1} failed ({type(exc).__name__})", flush=True)
                    time.sleep(2 * (attempt + 1))
            if result is None:
                result = {
                    "id": event_id,
                    "depth_type": "fetch_failed",
                    "origin_source": "",
                    "depth_km": "",
                    "origin_product_id": "",
                    "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                }
            writer.writerow(result)
            output.flush()
            if index % 50 == 0:
                print(f"{index}/{len(rows)}", flush=True)
            time.sleep(pause_seconds)
    return len(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_source = SNAPSHOT_V2 if os.path.exists(SNAPSHOT_V2) else SOURCE_LEGACY
    parser.add_argument("--source", default=default_source)
    parser.add_argument("--output", default=OUTPUT_DEFAULT)
    parser.add_argument("--pause-seconds", type=float, default=0.2)
    args = parser.parse_args(argv)
    count = fetch_flags(args.source, args.output, pause_seconds=args.pause_seconds)
    print(f"processed {count} source rows -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
