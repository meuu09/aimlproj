"""Merge the USGS ComCat and ISC Bulletin catalogs into data/combined_usgs_isc_catalog.csv."""
import os

import numpy as np
import pandas as pd

from seis_utils import (BOX, COMBINED_CSV, DATA, DEPTH_TYPES, SOLVED_DEPTH, USGS_CSV, flag_duplicates,
                        match_events, norm_magtype)

MATCH_S, MATCH_KM = 20.0, 100.0
COLS = ["id", "time", "latitude", "longitude", "horizontalError", "depth", "depthError", "depth_type", "depth_source",
        "location_source", "origin_id", "mag", "magType", "magType_norm", "mag_source", "source", "isc_id", "place",
        "match_dt_s", "match_dist_km", "isc_mag", "isc_magType", "usgs_time", "usgs_latitude", "usgs_longitude",
        "usgs_depth", "usgs_depth_type", "possible_duplicate_of"]


def load_isc():
    """ISC events in the box; sentinel errors (<= 0 or >= 999 km, e.g. NDI's 1000 km) set to missing."""
    isc = pd.read_csv(os.path.join(DATA, "isc_bulletin_events.csv"), dtype={"isc_id": str, "origin_id": str})
    isc = isc[isc.prime_latitude.between(*BOX["lat"]) & isc.prime_longitude.between(*BOX["lon"])]
    isc = isc.reset_index(drop=True)
    for col in ("horizontalError", "depthError"):
        isc.loc[(isc[col] <= 0) | (isc[col] >= 999), col] = np.nan
    return isc


def usgs_row(u):
    """Row from the USGS origin."""
    return dict(id=u.id, time=u.time, latitude=u.latitude, longitude=u.longitude, horizontalError=u.horizontalError,
                depth=u.depth, depthError=u.depthError, depth_type=u.depth_type, depth_source="USGS",
                location_source="USGS", origin_id=u.id, mag=u.mag, magType=u.magType, mag_source="USGS",
                source="USGS", isc_id=None, place=u.place)


def isc_origin(k):
    """Location fields of the whole ISC origin chosen in fetch_isc_bulletin.summarise."""
    return dict(time=k.time, latitude=k.latitude, longitude=k.longitude, horizontalError=k.horizontalError,
                depth=k.depth, depthError=k.depthError, depth_type=k.depth_type,
                depth_source=f"ISC:{k.origin_author}", location_source=f"ISC:{k.origin_author}",
                origin_id=k.origin_id)


def main():
    """Match, merge with one whole origin per row, flag possible duplicates, write the CSV."""
    us = pd.read_csv(USGS_CSV)
    us = us.merge(pd.read_csv(DEPTH_TYPES)[["id", "depth_type"]], on="id", how="left")
    us["depth_type"] = us.depth_type.fillna("unknown")
    isc = load_isc()
    prime = isc.rename(columns={"time": "_t", "latitude": "_la", "longitude": "_lo"}).rename(
        columns={"prime_time": "time", "prime_latitude": "latitude", "prime_longitude": "longitude"})
    m = match_events(us, prime, MATCH_S, MATCH_KM)
    match_u = {r.ia: r for r in m.itertuples()}
    matched_isc = set(m.ib)

    rows = []
    for i, u in us.iterrows():
        r = usgs_row(u)
        if i in match_u:
            mm = match_u[i]
            k = isc.loc[mm.ib]
            r.update(isc_id=k.isc_id, source="USGS+ISC", match_dt_s=mm.dt_s, match_dist_km=mm.dist_km,
                     isc_mag=k.mag, isc_magType=k.magType)
            usgs_solved = u.depth_type in SOLVED_DEPTH and (u.depthError or 0) > 0
            if not usgs_solved and k.depth_type == "from location":
                r.update(usgs_time=u.time, usgs_latitude=u.latitude, usgs_longitude=u.longitude,
                         usgs_depth=u.depth, usgs_depth_type=u.depth_type)
                r.update(isc_origin(k))
        rows.append(r)
    for j, k in isc.iterrows():
        if j in matched_isc or pd.isna(k.mag):
            continue
        r = dict(id=f"isc{k.isc_id}", mag=k.mag, magType=k.magType, mag_source=f"ISC:{k.mag_author}",
                 source="ISC", isc_id=k.isc_id, place=k.region)
        r.update(isc_origin(k))
        rows.append(r)
    out = pd.DataFrame(rows, columns=COLS)
    out["magType_norm"] = norm_magtype(out.magType)
    out["t"] = pd.to_datetime(out.time, utc=True, format="ISO8601")
    out = out.sort_values("t", kind="stable").drop(columns="t").reset_index(drop=True)
    out["possible_duplicate_of"], dup_pairs = flag_duplicates(out)
    out.to_csv(COMBINED_CSV, index=False)

    solved = out.depth_type.isin(SOLVED_DEPTH)
    print(f"USGS {len(us)} | ISC {len(isc)} | matched {len(m)} | ISC-only added {(out.source == 'ISC').sum()} "
          f"| combined {len(out)}")
    print("depth types:", out.depth_type.value_counts().to_dict())
    print("location sources (solved):", out[solved].location_source.value_counts().to_dict())
    print("matched USGS events given the whole ISC origin:", int(out.usgs_depth_type.notna().sum()))
    print("magType_norm:", out.magType_norm.value_counts().to_dict())
    print(f"possible duplicates flagged: {out.possible_duplicate_of.notna().sum()} rows "
          f"({int(dup_pairs.eligible.sum())} eligible close pairs)")


if __name__ == "__main__":
    main()
