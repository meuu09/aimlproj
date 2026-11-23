"""Fetch the USGS ComCat origin 'depth-type' flag for every event in our CSV.

depth-type tells whether the depth was solved for ("from location") or fixed by an
analyst ("operator assigned"). Output: data/usgs_depth_type.csv (id, depth_type, ...).
Run once:  .venv/bin/python fetch_depth_type.py   (about 5-10 minutes, 1 request per event)
"""
import csv, json, os, time, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "data", "usgs_himachal_kashmir_M3plus_1970-2026.csv")
DST = os.path.join(HERE, "data", "usgs_depth_type.csv")
URL = "https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={}&format=geojson"

done = {}
if os.path.exists(DST):
    done = {r["id"]: r for r in csv.DictReader(open(DST))}
rows = list(csv.DictReader(open(SRC)))
with open(DST, "a", newline="") as f:
    w = csv.DictWriter(f, ["id", "depth_type", "origin_source", "depth_km"])
    if not done:
        w.writeheader()
    for i, r in enumerate(rows):
        if r["id"] in done:
            continue
        for attempt in range(3):
            try:
                d = json.load(urllib.request.urlopen(URL.format(r["id"]), timeout=30))
                o = d["properties"]["products"]["origin"][0]["properties"]
                dt = (o.get("depth-type") or "").lower().replace("_", " ")
                w.writerow(dict(id=r["id"], depth_type=dt or "unknown",
                                origin_source=o.get("origin-source", ""), depth_km=o.get("depth", "")))
                break
            except Exception as e:
                time.sleep(2)
        else:
            w.writerow(dict(id=r["id"], depth_type="fetch_failed", origin_source="", depth_km=""))
        f.flush()
        if i % 50 == 0:
            print(i, "/", len(rows), flush=True)
        time.sleep(0.2)
print("done ->", DST)
