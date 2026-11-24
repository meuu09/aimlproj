"""Download the ISC Bulletin (ISF format) for the study box and parse it into event and origin tables."""
import csv
import os
import re
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "data", "isc_raw")
OUT_EVENTS = os.path.join(HERE, "data", "isc_bulletin_events.csv")
OUT_ORIGINS = os.path.join(HERE, "data", "isc_bulletin_origins.csv")
BOX = dict(bot_lat=31, top_lat=35, left_lon=74, right_lon=79.5)
YEARS = range(1970, 2027)
DEPTH_PRIORITY = ["ISC-EHB", "ISC", "NEIC", "NDI"]
MAG_PRIORITY = [("Mw", None), ("MW", None), ("mb", "ISC"), ("mb", "NEIC"), ("mb", None),
                ("ML", "NDI"), ("ML", None), ("MS", None), ("Ms", None)]
URL = ("https://www.isc.ac.uk/cgi-bin/web-db-run?request=COMPREHENSIVE&out_format=ISF"
       "&searchshape=RECT&bot_lat={bot_lat}&top_lat={top_lat}&left_lon={left_lon}&right_lon={right_lon}"
       "&start_year={y}&start_month=1&start_day=01&start_time=00%3A00%3A00"
       "&end_year={y}&end_month=12&end_day=31&end_time=23%3A59%3A59"
       "&req_mag_type=Any&req_mag_agcy=Any&include_magnitudes=on&include_headers=on&include_comments=on")
ORIGIN_RE = re.compile(r"^\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}")


def fnum(s):
    """Float from a fixed-width field, or None."""
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return None


def download(year):
    """Download one year of the ISF bulletin into data/isc_raw/ (skips existing files)."""
    path = os.path.join(RAW, f"isc_isf_{year}.txt")
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    for attempt in range(4):
        try:
            txt = urllib.request.urlopen(URL.format(y=year, **BOX), timeout=600).read().decode("utf-8", "replace")
            if "STOP" in txt or "No events were found" in txt:
                open(path, "w", encoding="utf-8").write(txt)
                return path
        except Exception as e:
            print(f"  {year}: attempt {attempt + 1} failed ({type(e).__name__})", flush=True)
        time.sleep(15 * (attempt + 1))
    raise RuntimeError(f"ISC download failed for {year}")


def parse(path):
    """Yield one dict per event: {'id','region','origins':[...],'mags':[...]}."""
    ev, last_origin, mode = None, None, None
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line.startswith("Event "):
            if ev:
                yield ev
            parts = line.split(None, 2)
            ev = dict(id=parts[1], region=parts[2] if len(parts) > 2 else "", origins=[], mags=[])
            mode = None
            continue
        if ev is None:
            continue
        if line.lstrip().startswith("Date       Time"):
            mode = "origin"; continue
        if line.startswith("Magnitude"):
            mode = "mag"; continue
        if mode == "origin" and ORIGIN_RE.match(line):
            o = dict(time=line[0:10].replace("/", "-") + "T" + line[11:22].strip(),
                     lat=fnum(line[36:44]), lon=fnum(line[45:54]),
                     depth=fnum(line[71:76]), depth_fixed=line[76:77] == "f",
                     depth_err=fnum(line[78:82]), ndef=fnum(line[83:87]),
                     author=line[118:127].strip(), origid=line[128:136].strip(), prime=False,
                     ehb_quality="")
            ev["origins"].append(o); last_origin = o
            continue
        if mode == "origin" and "(#PRIME)" in line and last_origin is not None:
            last_origin["prime"] = True; continue
        if mode == "origin" and "DEPTH_QUALITY=" in line and last_origin is not None:
            last_origin["ehb_quality"] = line.split("DEPTH_QUALITY=")[1].strip(" )"); continue
        if mode == "mag" and line.strip() and not line.startswith(" "):
            tok = line.split()
            if len(tok) >= 4 and fnum(tok[1]) is not None:
                ev["mags"].append(dict(type=tok[0], mag=fnum(tok[1]), author=tok[-2]))
    if ev:
        yield ev


def summarise(ev):
    """One row per event: prime epicentre, best solved depth and preferred magnitude."""
    ors = [o for o in ev["origins"] if o["lat"] is not None]
    if not ors:
        return None
    prime = next((o for o in ors if o["prime"]), ors[-1])
    free = [o for o in ors if o["depth"] is not None and not o["depth_fixed"] and (o["depth_err"] or 0) > 0]
    rank = lambda o: DEPTH_PRIORITY.index(o["author"]) if o["author"] in DEPTH_PRIORITY else len(DEPTH_PRIORITY)
    best = min(free, key=rank) if free else None
    if best:
        depth, derr, dtype, dauth = best["depth"], best["depth_err"], "from location", best["author"]
    else:
        depth, derr, dauth = prime["depth"], None, prime["author"]
        dtype = "operator assigned" if prime["depth_fixed"] else "unknown"
    mag = None
    for mt, ma in MAG_PRIORITY:
        mag = next((m for m in ev["mags"] if m["type"] == mt and (ma is None or m["author"] == ma)), None)
        if mag:
            break
    if mag is None and ev["mags"]:
        mag = ev["mags"][0]
    return dict(isc_id=ev["id"], time=prime["time"], latitude=prime["lat"], longitude=prime["lon"],
                prime_author=prime["author"], depth=depth, depthError=derr, depth_type=dtype,
                depth_author=dauth, mag=mag["mag"] if mag else None, magType=mag["type"] if mag else None,
                mag_author=mag["author"] if mag else None, n_origins=len(ors),
                agencies=";".join(sorted({o["author"] for o in ors})), region=ev["region"])


if __name__ == "__main__":
    os.makedirs(RAW, exist_ok=True)
    years = [int(sys.argv[1])] if len(sys.argv) > 1 else YEARS
    for y in years:
        download(y)
        print(f"{y}: downloaded", flush=True)
        time.sleep(3)
    events, origins = [], []
    for y in YEARS:
        p = os.path.join(RAW, f"isc_isf_{y}.txt")
        if not os.path.exists(p):
            continue
        for ev in parse(p):
            row = summarise(ev)
            if row:
                events.append(row)
            for o in ev["origins"]:
                origins.append(dict(isc_id=ev["id"], **o))
    for path, rows in ((OUT_EVENTS, events), (OUT_ORIGINS, origins)):
        if rows:
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"events={len(events)} origins={len(origins)} -> {OUT_EVENTS}")
