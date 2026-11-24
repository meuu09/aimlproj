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
Magnitudes are mixed: USGS mostly mb; ISC-only events mostly ML/MD from NDI.
Magnitude types are kept as reported (magType, mag_source); no conversion is made.

Extra columns (diagnostics, nothing is removed):
  horizontalError             USGS epicentre error (km) for USGS rows; empty for ISC-only rows
  match_dt_s, match_dist_km   time and distance between the matched USGS and ISC origins
  isc_mag, isc_magType        the ISC-chosen magnitude of a matched event, for comparison
  possible_duplicate_of       id of the row this one probably duplicates (same earthquake
                              listed twice, e.g. an NDI report the ISC did not associate):
                              within 60 s and 50 km, at least one of the two ISC-only.
                              See seis_utils.flag_duplicates. Analysis scripts drop these
                              rows only with --drop-duplicates.

Run after fetch_depth_type.py and fetch_isc_bulletin.py:
  .venv/bin/python build_combined_catalog.py
"""
import os
import numpy as np
import pandas as pd

from seis_utils import flag_duplicates, match_events

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
MATCH_S, MATCH_KM = 20.0, 100.0
TRUSTED_DEPTH = ("ISC-EHB", "ISC", "NEIC", "NDI")
SOLVED = ("from location", "constrained by depth phases")
COLS = ["id", "time", "latitude", "longitude", "horizontalError", "depth", "depthError", "depth_type", "depth_source",
        "mag", "magType", "mag_source", "source", "isc_id", "place",
        "match_dt_s", "match_dist_km", "isc_mag", "isc_magType", "possible_duplicate_of"]

us = pd.read_csv(os.path.join(D, "usgs_himachal_kashmir_M3plus_1970-2026.csv"))
us = us.merge(pd.read_csv(os.path.join(D, "usgs_depth_type.csv"))[["id", "depth_type"]], on="id", how="left")
us["depth_type"] = us.depth_type.fillna("unknown")
isc = pd.read_csv(os.path.join(D, "isc_bulletin_events.csv"), dtype={"isc_id": str})
isc = isc[isc.latitude.between(31, 35) & isc.longitude.between(74, 79.5)].reset_index(drop=True)
untrusted = (isc.depth_type == "from location") & ~isc.depth_author.isin(TRUSTED_DEPTH)
isc.loc[untrusted, ["depthError"]] = np.nan
isc.loc[untrusted, "depth_type"] = "unknown"

# greedy closest-first one-to-one matching (seis_utils.match_events)
m = match_events(us, isc, MATCH_S, MATCH_KM)
match_u = dict(zip(m.ia, m.ib))
match_i = dict(zip(m.ib, m.ia))
diag_u = {r.ia: (r.dt_s, r.dist_km) for r in m.itertuples()}

rows = []
for i, u in us.iterrows():
    r = dict(id=u.id, time=u.time, latitude=u.latitude, longitude=u.longitude,
             horizontalError=u.horizontalError, depth=u.depth,
             depthError=u.depthError, depth_type=u.depth_type, depth_source="USGS", mag=u.mag,
             magType=u.magType, mag_source="USGS", source="USGS", isc_id=None, place=u.place)
    if i in match_u:
        k = isc.loc[match_u[i]]
        r["isc_id"], r["source"] = k.isc_id, "USGS+ISC"
        r["match_dt_s"], r["match_dist_km"] = diag_u[i]
        r["isc_mag"], r["isc_magType"] = k.mag, k.magType
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
out = out.sort_values("t", kind="stable").drop(columns="t").reset_index(drop=True)
out["possible_duplicate_of"], dup_pairs = flag_duplicates(out)
out.to_csv(os.path.join(D, "combined_usgs_isc_catalog.csv"), index=False)

solved = out.depth_type.isin(SOLVED)
print(f"USGS {len(us)} | ISC {len(isc)} | matched {len(match_u)} | ISC-only added {(out.source == 'ISC').sum()} "
      f"| combined {len(out)}")
print("depth types:", out.depth_type.value_counts().to_dict())
print("depth sources (solved):", out[solved].depth_source.value_counts().to_dict())
print("USGS events whose depth was upgraded from ISC:", ((out.source == "USGS+ISC") & out.depth_source.str.startswith("ISC")).sum())
print("magTypes:", out.magType.value_counts().head(8).to_dict())
print(f"possible duplicates flagged: {out.possible_duplicate_of.notna().sum()} rows "
      f"({dup_pairs.flagged.sum()} flagged close pairs, {(~dup_pairs.flagged).sum()} USGS-USGS close pairs not flagged)")
