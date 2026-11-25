"""Download Global CMT solutions for the study box, derive P/T/B axes and faulting style, match to our catalog."""
import os
import re
import sys
import urllib.request

import numpy as np
import pandas as pd

from seis_utils import BOX, COMBINED_CSV, DATA, HERE, depth_class, match_events

RAW = os.path.join(DATA, "gcmt_box_cmtsolution.txt")
OUT = os.path.join(HERE, "results", "geology")
URL = ("https://www.globalcmt.org/cgi-bin/globalcmt-cgi-bin/CMT5/form?itype=ymd&yr=1976&mo=1&day=1&otype=ymd"
       "&oyr=2026&omo=11&oday=15&jyr=1976&jday=1&ojyr=1976&ojday=1&nday=1&lmw=0&umw=10&lms=0&ums=10&lmb=0&umb=10"
       f"&llat={BOX['lat'][0]}&ulat={BOX['lat'][1]}&llon={BOX['lon'][0]}&ulon={BOX['lon'][1]}&lhd=0&uhd=1000"
       "&lts=-9999&uts=9999&lpe1=0&upe1=90&lpe2=0&upe2=90&list=4")


def download():
    """Save the CMTSOLUTION text (HTML tags removed) unless it already exists."""
    if os.path.exists(RAW):
        return
    html = urllib.request.urlopen(URL, timeout=120).read().decode("utf-8", "replace")
    open(RAW, "w", encoding="utf-8").write(re.sub(r"<[^>]*>", "", html))


def parse(text):
    """One dict per CMTSOLUTION block: PDE hypocentre line plus centroid and moment tensor (dyne-cm)."""
    head = re.compile(r"^\s*([A-Z]{3,4})\s*(\d{4})\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+(-?[\d.]+)"
                      r"\s+(-?[\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)")
    rows, cur = [], None
    for line in text.splitlines():
        h = head.match(line)
        if h:
            p = h.groups()
            cur = dict(pde_time=pd.Timestamp(int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]))
                       + pd.Timedelta(seconds=float(p[6])), pde_latitude=float(p[7]), pde_longitude=float(p[8]),
                       pde_depth=float(p[9]), pde_mb=float(p[10]), pde_ms=float(p[11]), hypocentre_agency=p[0])
            rows.append(cur)
        elif cur is not None and ":" in line:
            k, v = [s.strip() for s in line.split(":", 1)]
            cur[k.replace(" ", "_")] = v
    return rows


def axes(m):
    """Plunge and azimuth of T, B, P (largest, middle, smallest eigenvalue) from r/t/p components."""
    mrr, mtt, mpp, mrt, mrp, mtp = m
    ned = np.array([[mtt, -mtp, mrt], [-mtp, mpp, -mrp], [mrt, -mrp, mrr]])
    ned[2, :] *= -1
    ned[:, 2] *= -1
    w, v = np.linalg.eigh(ned)
    out = {}
    for name, k in (("T", 2), ("B", 1), ("P", 0)):
        n, e, d = v[:, k] if v[2, k] >= 0 else -v[:, k]
        out[f"{name}_plunge"] = float(np.degrees(np.arcsin(np.clip(d, -1, 1))))
        out[f"{name}_azimuth"] = float(np.degrees(np.arctan2(e, n)) % 360)
    return out


def regime(p, b, t):
    """World Stress Map regime from axis plunges (Zoback 1992)."""
    if p >= 52 and t <= 35:
        return "normal"
    if 40 <= p < 52 and t <= 20:
        return "normal with strike-slip"
    if (p < 40 and b >= 45 and t <= 20) or (p <= 20 and b >= 45 and t < 40):
        return "strike-slip"
    if p <= 20 and 40 <= t < 52:
        return "thrust with strike-slip"
    if p <= 35 and t >= 52:
        return "thrust"
    return "unknown (oblique)"


def main():
    """Download, parse, classify, match with the combined catalog and write focal_mechanisms.csv."""
    os.makedirs(OUT, exist_ok=True)
    if "--parse-only" not in sys.argv:
        download()
    rows = parse(open(RAW, encoding="utf-8").read())
    out = []
    for r in rows:
        m = [float(r[k]) for k in ("Mrr", "Mtt", "Mpp", "Mrt", "Mrp", "Mtp")]
        m0 = np.sqrt(0.5 * sum(x * x for x in m[:3]) + sum(x * x for x in m[3:]))
        a = axes(m)
        out.append(dict(event_name=r["event_name"], pde_time=r["pde_time"], pde_latitude=r["pde_latitude"],
                        pde_longitude=r["pde_longitude"], pde_depth=r["pde_depth"],
                        centroid_latitude=float(r["latitude"]), centroid_longitude=float(r["longitude"]),
                        centroid_depth=float(r["depth"]), Mw=2 / 3 * np.log10(m0) - 10.7, **a,
                        regime=regime(a["P_plunge"], a["B_plunge"], a["T_plunge"]), source="Global CMT",
                        pde_in_box=bool(BOX["lat"][0] <= r["pde_latitude"] <= BOX["lat"][1]
                                        and BOX["lon"][0] <= r["pde_longitude"] <= BOX["lon"][1])))
    fm = pd.DataFrame(out)
    co = pd.read_csv(COMBINED_CSV)
    co["depth_class"] = depth_class(co.depth, co.depthError)
    mm = match_events(fm.rename(columns={"pde_time": "time", "pde_latitude": "latitude", "pde_longitude": "longitude"}),
                      co, 30, 100)
    fm["catalog_id"], fm["catalog_depth"], fm["catalog_depthError"], fm["catalog_depth_type"] = None, np.nan, np.nan, None
    fm["catalog_depth_class"] = None
    for ia, ib in zip(mm.ia, mm.ib):
        c = co.iloc[ib]
        fm.loc[ia, ["catalog_id", "catalog_depth", "catalog_depthError", "catalog_depth_type", "catalog_depth_class"]] = [
            c.id, c.depth, c.depthError, c.depth_type, c.depth_class]
    fm.to_csv(os.path.join(OUT, "focal_mechanisms.csv"), index=False)
    pd.set_option("display.width", 220)
    print(fm[["pde_time", "centroid_depth", "Mw", "P_plunge", "B_plunge", "T_plunge", "regime", "catalog_depth",
              "catalog_depthError", "catalog_depth_class"]].to_string(index=False, float_format="%.1f"))
    print(fm.regime.value_counts().to_string())


if __name__ == "__main__":
    main()
