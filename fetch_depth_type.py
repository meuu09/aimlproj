"""Fetch the USGS origin depth-type flag for every event into data/usgs_depth_type.csv."""
import csv
import json
import os
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "data", "usgs_himachal_kashmir_M3plus_1970-2026.csv")
DST = os.path.join(HERE, "data", "usgs_depth_type.csv")
LOG = os.path.join(HERE, "data", "usgs_depth_type_failures.log")
URL = "https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={}&format=geojson"
FIELDS = ["id", "depth_type", "origin_source", "depth_km"]


def fetch(event_id, attempts=4):
    """Depth-type row for one event, retrying with exponential back-off; None after the last failure."""
    for k in range(attempts):
        try:
            d = json.load(urllib.request.urlopen(URL.format(event_id), timeout=30))
            o = d["properties"]["products"]["origin"][0]["properties"]
            dt = (o.get("depth-type") or "").lower().replace("_", " ")
            return dict(id=event_id, depth_type=dt or "unknown", origin_source=o.get("origin-source", ""),
                        depth_km=o.get("depth", ""))
        except Exception as e:
            with open(LOG, "a", encoding="utf-8") as f:
                f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%S')} {event_id} attempt {k + 1}: {type(e).__name__}: {e}\n")
            time.sleep(2 ** k)
    return None


def save(rows):
    """Rewrite the output file (one row per event)."""
    with open(DST, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, FIELDS)
        w.writeheader()
        w.writerows(rows.values())


def main():
    """Fetch every event that is missing or failed before; successful rows are not fetched again."""
    done = {}
    if os.path.exists(DST):
        done = {r["id"]: r for r in csv.DictReader(open(DST, encoding="utf-8"))}
    ids = [r["id"] for r in csv.DictReader(open(SRC, encoding="utf-8"))]
    todo = [i for i in ids if done.get(i, {}).get("depth_type", "fetch_failed") == "fetch_failed"]
    print(f"{len(ids)} events, {len(todo)} to fetch")
    for n, event_id in enumerate(todo, 1):
        done[event_id] = fetch(event_id) or dict(id=event_id, depth_type="fetch_failed", origin_source="",
                                                 depth_km="")
        if n % 50 == 0:
            save(done)
            print(n, "/", len(todo), flush=True)
        time.sleep(0.2)
    save({i: done[i] for i in ids if i in done})
    failed = sum(r["depth_type"] == "fetch_failed" for r in done.values())
    print(f"done -> {DST}; failed: {failed}" + (f" (see {LOG})" if failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
