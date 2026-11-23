"""
Merge the USGS ComCat download with the ISC Bulletin into one catalog of real events
(no synthetic data): data/combined_usgs_isc_catalog.csv

Matching: an ISC event and a USGS event are the same earthquake if their origin
times differ by <= MATCH_S seconds and epicentres by <= MATCH_KM km (closest match,
each event used once).

Per merged event:
  * matched     USGS row (id, time, place, magnitude) kept. Depth: USGS depth if USGS
                solved it (depth-type from location / depth phases, error > 0);
                otherwise the ISC-chosen depth if that came from a trusted agency
                with an error; otherwise the USGS (fixed/unflagged) depth.
  * ISC only    added from the ISC Bulletin (mostly small events located by India's
                National Centre for Seismology, NDI).
  * USGS only   kept as is.
Depths from agencies outside TRUSTED_DEPTH are treated as unknown.
Magnitudes are mixed: USGS mostly mb; ISC-only events mostly ML from NDI.

Run after fetch_depth_type.py and fetch_isc_bulletin.py:
  .venv/bin/python build_combined_catalog.py
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
MATCH_S, MATCH_KM = 20.0, 100.0
TRUSTED_DEPTH = ("ISC-EHB", "ISC", "NEIC", "NDI")
SOLVED = ("from location", "constrained by depth phases")
COLS = ["id", "time", "latitude", "longitude", "depth", "depthError", "depth_type", "depth_source",
        "mag", "magType", "mag_source", "source", "isc_id", "place"]

us = pd.read_csv(os.path.join(D, "usgs_himachal_kashmir_M3plus_1970-2026.csv"))
us = us.merge(pd.read_csv(os.path.join(D, "usgs_depth_type.csv"))[["id", "depth_type"]], on="id", how="left")
us["depth_type"] = us.depth_type.fillna("unknown")
isc = pd.read_csv(os.path.join(D, "isc_bulletin_events.csv"), dtype={"isc_id": str})
isc = isc[isc.latitude.between(31, 35) & isc.longitude.between(74, 79.5)].reset_index(drop=True)
untrusted = (isc.depth_type == "from location") & ~isc.depth_author.isin(TRUSTED_DEPTH)
isc.loc[untrusted, ["depthError"]] = np.nan
isc.loc[untrusted, "depth_type"] = "unknown"

tu = pd.to_datetime(us.time, utc=True, format="ISO8601").values.astype("datetime64[ms]").astype(float) / 1000
ti = pd.to_datetime(isc.time, utc=True, format="ISO8601").values.astype("datetime64[ms]").astype(float) / 1000


def dist_km(lat1, lon1, lat2, lon2):
    p = np.radians
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


# candidate pairs, then greedy closest-first one-to-one matching
pairs = []
for i in range(len(us)):
    near = np.where(np.abs(ti - tu[i]) <= MATCH_S)[0]
    for j in near:
        d = dist_km(us.latitude[i], us.longitude[i], isc.latitude[j], isc.longitude[j])
        if d <= MATCH_KM:
            pairs.append((abs(ti[j] - tu[i]) / MATCH_S + d / MATCH_KM, i, j))
match_u, match_i = {}, {}
for _, i, j in sorted(pairs):
    if i not in match_u and j not in match_i:
        match_u[i], match_i[j] = j, i

rows = []
for i, u in us.iterrows():
    r = dict(id=u.id, time=u.time, latitude=u.latitude, longitude=u.longitude, depth=u.depth,
             depthError=u.depthError, depth_type=u.depth_type, depth_source="USGS", mag=u.mag,
             magType=u.magType, mag_source="USGS", source="USGS", isc_id=None, place=u.place)
    if i in match_u:
        k = isc.loc[match_u[i]]
        r["isc_id"], r["source"] = k.isc_id, "USGS+ISC"
        usgs_solved = u.depth_type in SOLVED and (u.depthError or 0) > 0
        if not usgs_solved and k.depth_type == "from location":
            r.update(depth=k.depth, depthError=k.depthError, depth_type="from location",
                     depth_source=f"ISC:{k.depth_author}")
    rows.append(r)
for j, k in isc.iterrows():
    if j in match_i or pd.isna(k.mag):
        continue
    rows.append(dict(id=f"isc{k.isc_id}", time=k.time, latitude=k.latitude, longitude=k.longitude,
                     depth=k.depth, depthError=k.depthError, depth_type=k.depth_type,
                     depth_source=f"ISC:{k.depth_author}", mag=k.mag, magType=k.magType,
                     mag_source=f"ISC:{k.mag_author}", source="ISC", isc_id=k.isc_id, place=k.region))
out = pd.DataFrame(rows, columns=COLS)
out["t"] = pd.to_datetime(out.time, utc=True, format="ISO8601")
out = out.sort_values("t").drop(columns="t")
out.to_csv(os.path.join(D, "combined_usgs_isc_catalog.csv"), index=False)

solved = out.depth_type.isin(SOLVED)
print(f"USGS {len(us)} | ISC {len(isc)} | matched {len(match_u)} | ISC-only added {(out.source == 'ISC').sum()} "
      f"| combined {len(out)}")
print("depth types:", out.depth_type.value_counts().to_dict())
print("depth sources (solved):", out[solved].depth_source.value_counts().to_dict())
print("USGS events whose depth was upgraded from ISC:", ((out.source == "USGS+ISC") & out.depth_source.str.startswith("ISC")).sum())
print("magTypes:", out.magType.value_counts().head(8).to_dict())
