"""Download a dated USGS ComCat CSV snapshot to a new, versioned file.

The repository's existing snapshot remains untouched. The default output is
``data/usgs_himachal_kashmir_M3plus_1970-2026_snapshot_v2.csv`` and the exact
query is written next to it as JSON.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import urllib.parse
import urllib.request
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTPUT = os.path.join(HERE, "data", "usgs_himachal_kashmir_M3plus_1970-2026_snapshot_v2.csv")
API = "https://earthquake.usgs.gov/fdsnws/event/1/query"
DEFAULT_QUERY = {
    "format": "csv",
    "starttime": "1970-03-05T00:00:00",
    "endtime": "2026-11-14T23:59:59",
    "minlatitude": "31.0",
    "maxlatitude": "35.0",
    "minlongitude": "74.0",
    "maxlongitude": "79.5",
    "minmagnitude": "3.0",
    "orderby": "time",
}


def query_url(parameters: dict[str, str]) -> str:
    return API + "?" + urllib.parse.urlencode(parameters)


def fetch(output_path: str, parameters: dict[str, str], *, force: bool = False) -> int:
    if os.path.exists(output_path) and not force:
        raise FileExistsError(f"Refusing to overwrite existing file: {output_path}; pass --force to replace it")
    request_url = query_url(parameters)
    request = urllib.request.Request(request_url, headers={"User-Agent": "himalayan-seismicity-research/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:
        content = response.read()
    reader = csv.DictReader(io.StringIO(content.decode("utf-8-sig")))
    if not reader.fieldnames or not {"id", "time", "latitude", "longitude", "depth", "mag", "magType"}.issubset(reader.fieldnames):
        raise ValueError("USGS response did not contain the expected ComCat CSV columns")
    count = sum(1 for _ in reader)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "wb") as output:
        output.write(content)
    metadata = {
        "source": "USGS FDSN Event Web Service / ComCat",
        "url": request_url,
        "parameters": parameters,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "rows": count,
    }
    with open(output_path + ".query.json", "w", encoding="utf-8") as sidecar:
        json.dump(metadata, sidecar, indent=2)
        sidecar.write("\n")
    print(f"Downloaded {count} ComCat rows to {output_path}")
    print(f"Query metadata: {output_path}.query.json")
    return count


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--starttime", default=DEFAULT_QUERY["starttime"])
    parser.add_argument("--endtime", default=DEFAULT_QUERY["endtime"])
    parser.add_argument("--minmagnitude", default=DEFAULT_QUERY["minmagnitude"])
    parser.add_argument("--minlatitude", default=DEFAULT_QUERY["minlatitude"])
    parser.add_argument("--maxlatitude", default=DEFAULT_QUERY["maxlatitude"])
    parser.add_argument("--minlongitude", default=DEFAULT_QUERY["minlongitude"])
    parser.add_argument("--maxlongitude", default=DEFAULT_QUERY["maxlongitude"])
    parser.add_argument("--force", action="store_true", help="allow replacing the specified output file")
    args = parser.parse_args(argv)
    parameters = {
        "format": "csv",
        "starttime": args.starttime,
        "endtime": args.endtime,
        "minlatitude": args.minlatitude,
        "maxlatitude": args.maxlatitude,
        "minlongitude": args.minlongitude,
        "maxlongitude": args.maxlongitude,
        "minmagnitude": args.minmagnitude,
        "orderby": "time",
    }
    fetch(args.output, parameters, force=args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
