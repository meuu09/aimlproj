"""Shared helpers: catalog loading, depth classes, Mc, b-values, event matching and duplicate flags."""
import argparse
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
USGS_CSV = os.path.join(DATA, "usgs_himachal_kashmir_M3plus_1970-2026.csv")
DEPTH_TYPES = os.path.join(DATA, "usgs_depth_type.csv")
COMBINED_CSV = os.path.join(DATA, "combined_usgs_isc_catalog.csv")

BOX = dict(lat=(31.0, 35.0), lon=(74.0, 79.5))
CATALOG_START = pd.Timestamp("1970-01-01")
CATALOG_END = pd.Timestamp("2026-11-15")
SOLVED_DEPTH = ("from location", "constrained by depth phases")
DEEP_MIN, DEEP_MAX = 40.0, 150.0
EARTH_R = 6371.0

CLUSTER_POLY = [(75.50, 33.60), (76.10, 33.60), (76.50, 32.95),
                (76.50, 32.75), (75.90, 32.75), (75.50, 33.20)]

DUP_MAX_S, DUP_MAX_KM = 60.0, 50.0

LONG_CLASS = {
    "intermediate": "intermediate (40 <= depth - error, depth + error <= 150 km)",
    "shallow": "shallow (depth + error < 40 km)",
    "uncertain": "uncertain (error bar crosses 40 or 150 km, or no error)",
    "deep": "deep (depth - error > 150 km)",
}
MAGTYPE_NORM = {"mb": "mb", "ml": "ML", "md": "MD", "ms": "Ms", "mw": "Mw", "mww": "Mw", "mwc": "Mw",
                "mwb": "Mw", "mwr": "Mw", "mi": "Mi", "mpv": "mpv"}


def parse_args(description="", argv=None):
    """Options shared by the analysis scripts him_ides_himachal, kishtwar_cluster and depth_clusters."""
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--catalog", choices=["usgs", "combined"], default="usgs",
                   help="usgs = USGS ComCat only; combined = USGS + ISC Bulletin")
    p.add_argument("--include-fixed", action="store_true",
                   help="keep fixed (analyst-assigned) and unflagged depths")
    p.add_argument("--drop-duplicates", action="store_true",
                   help="combined catalog only: drop rows flagged in possible_duplicate_of")
    a = p.parse_args(argv)
    if a.drop_duplicates and a.catalog != "combined":
        p.error("--drop-duplicates only applies to --catalog combined")
    return a


def results_base(args):
    """Folder that holds this catalog's case folders."""
    if args.catalog == "usgs":
        return os.path.join(HERE, "results")
    name = "combined_usgs_isc_nodup" if args.drop_duplicates else "combined_usgs_isc"
    return os.path.join(HERE, "results", name)


def case_name(args):
    """Name of the depth case folder."""
    return "including_fixed" if args.include_fixed else "excluding_fixed"


def case_label(args):
    """Catalog and depth case, as stamped on the figures."""
    cat = "USGS" if args.catalog == "usgs" else "USGS + ISC"
    if args.drop_duplicates:
        cat += " (possible duplicates removed)"
    case = ("Case: including fixed/unflagged depths" if args.include_fixed
            else "Case: excluding fixed depths (solved only)")
    return f"{cat} | {case}"


def load_catalog(catalog="usgs", include_fixed=False, start_year=None, drop_duplicates=False,
                 box=BOX, verbose=True):
    """Load a catalog, cut it to the box, start year and depth case, and add depth classes.
    Returns (df, info); info holds the event count left after each filtering step."""
    if catalog == "usgs":
        df = pd.read_csv(USGS_CSV)
        df = df.merge(pd.read_csv(DEPTH_TYPES)[["id", "depth_type"]], on="id", how="left")
        df["source"] = "USGS"
        df["depth_source"] = "USGS"
        df["possible_duplicate_of"] = np.nan
    elif catalog == "combined":
        df = pd.read_csv(COMBINED_CSV, dtype={"isc_id": str})
        if "possible_duplicate_of" not in df:
            raise RuntimeError("combined catalog has no possible_duplicate_of column: "
                               "rerun python -m catalog.build_combined_catalog")
    else:
        raise ValueError(f"unknown catalog {catalog!r}")
    if "horizontalError" not in df:
        df["horizontalError"] = np.nan
    df["depth_type"] = df.depth_type.fillna("unknown")
    df["time"] = pd.to_datetime(df.time, utc=True, format="ISO8601").dt.tz_localize(None)
    if "magType_norm" not in df:
        df["magType_norm"] = norm_magtype(df.magType)
    info = {"catalog": catalog, "in_file": len(df)}
    df = df[df.latitude.between(*box["lat"]) & df.longitude.between(*box["lon"]) & df.mag.notna()
            & (df.time >= CATALOG_START) & (df.time < CATALOG_END)]
    info["in_box_with_mag"] = len(df)
    if start_year is not None:
        df = df[df.time.dt.year >= start_year]
    info["after_start_year"] = len(df)
    if drop_duplicates:
        df = df[df.possible_duplicate_of.isna()]
    info["after_duplicates"] = len(df)
    if not include_fixed:
        df = df[df.depth_type.isin(SOLVED_DEPTH)]
    info["after_depth_filter"] = len(df)
    df = df.copy()
    df["depth_class"] = depth_class(df.depth, df.depthError)
    df = df.sort_values("time").reset_index(drop=True)
    if verbose:
        print(f"[load] {catalog}: file {info['in_file']} -> box {info['in_box_with_mag']} -> "
              f"{'all years' if start_year is None else f'year>={start_year}'} {info['after_start_year']} -> "
              f"dedup {info['after_duplicates']} -> "
              f"depth filter {info['after_depth_filter']}")
    return df, info


def depth_class(depth, depth_error, deep_min=DEEP_MIN, deep_max=DEEP_MAX):
    """Depth class from the whole error bar: intermediate if [depth - error, depth + error] lies in
    [40, 150] km, shallow if above 40 km, deep if below 150 km; otherwise (or no error) uncertain."""
    z = np.asarray(depth, dtype=float)
    e = np.asarray(depth_error, dtype=float)
    with np.errstate(invalid="ignore"):
        inter = (z - e >= deep_min) & (z + e <= deep_max)
        shal = (z + e) < deep_min
        deep = (z - e) > deep_max
    return np.select([inter, shal, deep], ["intermediate", "shallow", "deep"], default="uncertain")


def max_curv(mags, step=0.1):
    """Mc by maximum curvature (Wiemer & Wyss 2000): the most frequent magnitude bin.
    Bins are floor-rounded multiples of `step` so float noise cannot shift a magnitude."""
    m = np.asarray(mags, dtype=float)
    m = m[np.isfinite(m)]
    if m.size == 0:
        return np.nan
    k = np.floor(m / step + 1e-6).astype(int)
    vals, counts = np.unique(k, return_counts=True)
    return round(float(vals[np.argmax(counts)] * step), 2)


def max_curv_eqcorrscan(mags, bin_size=0.5):
    """The authors' Mc: eqcorrscan 0.4.2 calc_max_curv (peak |second derivative| of the
    0.5-bin histogram). The paper adds +0.4 to it."""
    min_bin, max_bin = int(min(mags)), int(max(mags) + 1)
    bins = np.arange(min_bin, max_bin + bin_size, bin_size)
    df, bins = np.histogram(mags, bins)
    curvature = np.diff(np.diff(df) / bin_size) / bin_size
    return float(bins[np.argmax(np.abs(curvature))] + bin_size)


def b_lsq(mags, mc, min_steps=3):
    """The repo's b-value: line fit to log10(cumulative number) for M >= mc.
    Returns (b, steps, cum, fit_fn); NaN and None with fewer than `min_steps` magnitudes."""
    m = np.asarray(mags, dtype=float)
    m = np.sort(m[m >= mc - 1e-9])
    steps = np.unique(m)
    cum = np.array([(m >= s - 1e-9).sum() for s in steps])
    if len(steps) < min_steps:
        return np.nan, steps, cum, None
    fit = np.polyfit(steps, np.log10(cum), 1)
    return -fit[0], steps, cum, np.poly1d(fit)


def b_mle(mags, mc, dm=0.1, min_n=5):
    """Aki (1965) / Utsu maximum-likelihood b-value with binning correction.
    Returns (b, Aki sigma b/sqrt(N), N); NaN below min_n events."""
    m = np.asarray(mags, dtype=float)
    m = m[m >= mc - 1e-9]
    if len(m) < min_n:
        return np.nan, np.nan, len(m)
    b = np.log10(np.e) / (m.mean() - (mc - dm / 2))
    return b, b / np.sqrt(len(m)), len(m)


def b_shi_bolt_sigma(mags, mc, dm=0.1):
    """Shi & Bolt (1982) 1-sigma of the MLE b-value, from the spread of the magnitudes."""
    m = np.asarray(mags, dtype=float)
    m = m[m >= mc - 1e-9]
    b, _, n = b_mle(m, mc, dm)
    if np.isnan(b):
        return np.nan
    return 2.30 * b ** 2 * np.sqrt(((m - m.mean()) ** 2).sum() / (n * (n - 1)))


def b_bootstrap(mags, mc, dm=0.1, n_boot=2000, seed=0, ci=95):
    """Bootstrap the MLE b-value: resample the events above mc with replacement.
    Returns (median, lower, upper) of the central `ci` percent; NaNs if too few events."""
    m = np.asarray(mags, dtype=float)
    m = m[m >= mc - 1e-9]
    if len(m) < 5:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    samples = rng.choice(m, size=(n_boot, len(m)), replace=True)
    spread = samples.mean(axis=1) - (mc - dm / 2)
    with np.errstate(divide="ignore"):
        bs = np.log10(np.e) / spread
    bs = bs[np.isfinite(bs)]
    lo, hi = np.percentile(bs, [(100 - ci) / 2, 100 - (100 - ci) / 2])
    return float(np.median(bs)), float(lo), float(hi)


def _eqcorrscan_or_nan(mags):
    """max_curv_eqcorrscan, or NaN when the magnitudes span too few bins for it."""
    try:
        return max_curv_eqcorrscan(mags)
    except ValueError:
        return np.nan


def mc_sliding(times, mags, window=50, step=10, n_boot=200, seed=0):
    """Mc in sliding windows of `window` consecutive events, moved by `step` events:
    max curvature with a bootstrap 16-84 % range, and the authors' Mc + 0.4."""
    order = np.argsort(np.asarray(times))
    t = pd.to_datetime(pd.Series(np.asarray(times)[order])).reset_index(drop=True)
    m = np.asarray(mags, dtype=float)[order]
    rng = np.random.default_rng(seed)
    rows = []
    for start in range(0, len(m) - window + 1, step):
        w = m[start:start + window]
        boot = [max_curv(rng.choice(w, len(w), replace=True)) for _ in range(n_boot)]
        rows.append(dict(t_start=t[start], t_end=t[start + window - 1],
                         t_mid=t[start] + (t[start + window - 1] - t[start]) / 2,
                         n=window, mc_maxc=max_curv(w), mc_maxc_p16=np.percentile(boot, 16),
                         mc_maxc_p84=np.percentile(boot, 84),
                         mc_authors=_eqcorrscan_or_nan(w) + 0.4))
    return pd.DataFrame(rows)


def interevent_days(times):
    """Inter-event times (days) of a time series, sorted first. Returns (times[1:], gaps)."""
    t = np.sort(pd.to_datetime(pd.Series(times)).values)
    return pd.to_datetime(t[1:]), np.diff(t).astype("timedelta64[s]").astype(float) / 86400.0


def cv(gaps):
    """Coefficient of variation of inter-event times (1 = Poisson, >1 = clustered)."""
    g = np.asarray(gaps, dtype=float)
    return g.std() / g.mean() if len(g) > 1 and g.mean() > 0 else np.nan


def norm_magtype(types):
    """Group equivalent magnitude labels (e.g. MW, mww, mwc -> Mw); unknown labels are kept as given."""
    s = pd.Series(types, dtype=object)
    return s.map(lambda t: MAGTYPE_NORM.get(str(t).lower(), t) if pd.notna(t) else "none")


def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance in km."""
    p = np.radians
    a = (np.sin(p(np.asarray(lat2) - lat1) / 2) ** 2
         + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(np.asarray(lon2) - lon1) / 2) ** 2)
    return EARTH_R * 2 * np.arcsin(np.sqrt(a))


def xy_km(lon, lat, lon0, lat0):
    """Flat-earth x/y in km around (lon0, lat0)."""
    return (lon - lon0) * 111.32 * np.cos(np.radians(lat0)), (lat - lat0) * 110.57


def _epoch_s(times):
    """Times as float seconds since 1970 (millisecond resolution)."""
    t = pd.to_datetime(pd.Series(times), utc=True, format="ISO8601").dt.tz_convert(None).to_numpy()
    return t.astype("datetime64[ms]").astype(float) / 1000


def match_events(a, b, max_s, max_km):
    """Greedy one-to-one matching of catalogs a and b within max_s seconds and max_km km,
    best score (dt/max_s + dist/max_km) first. Returns positions ia, ib with dt_s, dist_km."""
    ta, tb = _epoch_s(a.time), _epoch_s(b.time)
    order = np.argsort(tb)
    tb_sorted = tb[order]
    lat_a, lon_a = a.latitude.to_numpy(float), a.longitude.to_numpy(float)
    lat_b, lon_b = b.latitude.to_numpy(float), b.longitude.to_numpy(float)
    pairs = []
    for i in range(len(a)):
        lo = np.searchsorted(tb_sorted, ta[i] - max_s, side="left")
        hi = np.searchsorted(tb_sorted, ta[i] + max_s, side="right")
        for j in order[lo:hi]:
            d = float(haversine_km(lat_a[i], lon_a[i], lat_b[j], lon_b[j]))
            if d <= max_km:
                dt = abs(tb[j] - ta[i])
                pairs.append((dt / max_s + d / max_km, i, j, dt, d))
    used_a, used_b, out = set(), set(), []
    for _, i, j, dt, d in sorted(pairs):
        if i not in used_a and j not in used_b:
            used_a.add(i); used_b.add(j)
            out.append(dict(ia=i, ib=j, dt_s=dt, dist_km=d))
    return pd.DataFrame(out, columns=["ia", "ib", "dt_s", "dist_km"])


def close_pairs(df, max_s=DUP_MAX_S, max_km=DUP_MAX_KM):
    """All pairs of rows within max_s seconds and max_km km (positions in df)."""
    t = _epoch_s(df.time)
    order = np.argsort(t, kind="stable")
    lat, lon = df.latitude.to_numpy(float), df.longitude.to_numpy(float)
    out = []
    for n, i in enumerate(order):
        for j in order[n + 1:]:
            if t[j] - t[i] > max_s:
                break
            d = float(haversine_km(lat[i], lon[i], lat[j], lon[j]))
            if d <= max_km:
                out.append(dict(i=int(i), j=int(j), dt_s=float(t[j] - t[i]), dist_km=d))
    return pd.DataFrame(out, columns=["i", "j", "dt_s", "dist_km"])


def flag_duplicates(df, max_s=DUP_MAX_S, max_km=DUP_MAX_KM):
    """Flag ISC-only rows that probably repeat another row (within max_s and max_km). USGS rows are never
    flagged, so no chain can join two USGS events; ISC-only groups keep their best-quality row.
    Returns (possible_duplicate_of aligned to df.index, close pairs with 'eligible' and 'flagged')."""
    df = df.reset_index(drop=False).rename(columns={"index": "_orig_index"})
    pairs = close_pairs(df, max_s, max_km)
    isc = (df.source.astype(str) == "ISC").to_numpy()
    ii, jj = pairs.i.to_numpy(int), pairs.j.to_numpy(int)
    pairs["eligible"] = isc[ii] | isc[jj] if len(pairs) else pd.Series(dtype=bool)
    score = (pairs.dt_s / max_s + pairs.dist_km / max_km).to_numpy()
    dup = {}
    best_usgs = {}
    for i, j, s, ok in zip(ii, jj, score, pairs.eligible):
        for a, b in ((i, j), (j, i)):
            if ok and isc[a] and not isc[b] and s < best_usgs.get(a, (np.inf, None))[0]:
                best_usgs[a] = (s, b)
    for a, (_, b) in best_usgs.items():
        dup[a] = b
    parent = {k: k for k in np.flatnonzero(isc)}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in zip(ii, jj):
        if isc[i] and isc[j]:
            parent[find(j)] = find(i)
    groups = {}
    for k in parent:
        groups.setdefault(find(k), []).append(k)
    solved = df.depth_type.isin(SOLVED_DEPTH).to_numpy()
    has_err = (df.depthError.fillna(0) > 0).to_numpy()
    mag = df.mag.fillna(-9).to_numpy(float)
    tsec = _epoch_s(df.time)
    for members in groups.values():
        if len(members) < 2 or any(k in dup for k in members):
            continue
        best = min(members, key=lambda k: (not solved[k], not has_err[k], -mag[k], tsec[k]))
        for k in members:
            if k != best:
                dup[k] = best
    ids = df.id.to_numpy()
    dup_of = pd.Series(np.nan, index=df._orig_index, dtype=object)
    for k, b in dup.items():
        dup_of.iloc[k] = ids[b]
    pairs["flagged"] = [dup.get(i) == j or dup.get(j) == i for i, j in zip(ii, jj)]
    pairs.insert(0, "id_1", ids[ii] if len(pairs) else [])
    pairs.insert(1, "id_2", ids[jj] if len(pairs) else [])
    return dup_of, pairs.drop(columns=["i", "j"])


def b_reliability(n):
    """Label for a b-value from N events: Aki's sigma = b/sqrt(N) is ~14 % of b at N=50 and ~22 % at N=20."""
    return "reliable (N>=50)" if n >= 50 else ("rough (20<=N<50)" if n >= 20 else "not reported (N<20)")


def mc_gft(mags, dm=0.1, level=90.0):
    """Goodness-of-fit Mc (Wiemer & Wyss 2000): lowest Mc whose G-R fit explains >= `level` % of the
    observed non-cumulative counts above it (R = 100 - 100 sum|B-S| / sum B). NaN if never reached."""
    m = np.asarray(mags, dtype=float)
    m = m[np.isfinite(m)]
    if m.size < 20:
        return np.nan
    k = np.floor(m / dm + 1e-6).astype(int)
    for c in np.arange(k.min(), k.max() + 1):
        above = k[k >= c]
        if above.size < 20:
            break
        b, _, n = b_mle(above * dm, c * dm, dm)
        bins = np.arange(c, above.max() + 1)
        obs = np.array([(above == x).sum() for x in bins], dtype=float)
        cum = n * 10 ** (-b * (bins - c) * dm)
        syn = cum - np.append(cum[1:], 0.0)
        if 100 - 100 * np.abs(obs - syn).sum() / obs.sum() >= level:
            return round(float(c * dm), 2)
    return np.nan


def mc_mbs(mags, dm=0.1, width=0.5, min_n=50):
    """b-value-stability Mc (Cao & Gao 2002, as in Woessner & Wiemer 2005): lowest Mc where the mean b over
    [Mc, Mc + width] differs from b(Mc) by no more than its Shi-Bolt sigma. NaN if never reached."""
    m = np.asarray(mags, dtype=float)
    m = m[np.isfinite(m)]
    if m.size < min_n:
        return np.nan
    lo, hi = np.floor(m.min() / dm + 1e-6), np.floor(m.max() / dm + 1e-6)
    steps = int(round(width / dm))
    for c in np.arange(lo, hi + 1):
        mc = c * dm
        if (m >= mc - 1e-9).sum() < min_n:
            break
        b = b_mle(m, mc, dm)[0]
        avg = np.nanmean([b_mle(m, (c + i) * dm, dm)[0] for i in range(steps + 1)])
        if abs(avg - b) <= b_shi_bolt_sigma(m, mc, dm):
            return round(float(mc), 2)
    return np.nan


def poisson_ci(n, years, conf=0.95):
    """Exact (Garwood 1936) confidence interval of a Poisson rate n / years."""
    from scipy.stats import chi2
    a = 1 - conf
    lo = 0.0 if n == 0 else chi2.ppf(a / 2, 2 * n) / 2
    hi = chi2.ppf(1 - a / 2, 2 * n + 2) / 2
    return lo / years, hi / years


def cv_poisson_test(times, start, end, n_sim=2000, seed=0):
    """CV of inter-event times against a stationary Poisson process with the same number of events in the
    same window (uniform times). Returns (cv, 2.5 %, 97.5 % of simulated CV, two-sided p); NaN below 4 events."""
    t = np.sort(np.asarray(pd.to_datetime(pd.Series(times)).values.astype("datetime64[s]").astype(float)))
    if t.size < 4:
        return np.nan, np.nan, np.nan, np.nan
    obs = cv(np.diff(t))
    s0, s1 = [pd.Timestamp(x).to_datetime64().astype("datetime64[s]").astype(float) for x in (start, end)]
    rng = np.random.default_rng(seed)
    sims = np.sort(rng.uniform(s0, s1, size=(n_sim, t.size)), axis=1)
    g = np.diff(sims, axis=1)
    sim_cv = g.std(axis=1) / g.mean(axis=1)
    p = 2 * min((sim_cv >= obs).mean(), (sim_cv <= obs).mean())
    return float(obs), float(np.percentile(sim_cv, 2.5)), float(np.percentile(sim_cv, 97.5)), float(min(p, 1.0))


def dispersion_index(times, start, end):
    """Variance/mean of annual counts and its chi-square p-value against Poisson (overdispersion test)."""
    from scipy.stats import chi2
    years = np.arange(pd.Timestamp(start).year, pd.Timestamp(end).year + 1)
    counts = pd.Series(pd.to_datetime(pd.Series(times)).dt.year).value_counts().reindex(years, fill_value=0).values
    if len(counts) < 2 or counts.mean() == 0:
        return np.nan, np.nan
    d = counts.var(ddof=1) / counts.mean()
    return float(d), float(chi2.sf(d * (len(counts) - 1), len(counts) - 1))


def hypo_dist_km(lat1, lon1, z1, lat2, lon2, z2):
    """3-D hypocentral distance (km): great-circle horizontal distance combined with the depth difference."""
    return np.sqrt(haversine_km(lat1, lon1, lat2, lon2) ** 2 + (np.asarray(z2) - z1) ** 2)


def candidate_sequences(df, days=20.0, km=20.0, min_events=3):
    """Every group of an initiating event plus all later events within `days` and `km` (3-D) of it,
    as in Michailos et al. 2021 Sect. 3.2.3. All candidates are returned; overlaps are not removed."""
    d = df.sort_values("time").reset_index(drop=True)
    t = d.time.values
    out = []
    for i in range(len(d)):
        later = np.flatnonzero((t > t[i]) & (t <= t[i] + np.timedelta64(int(days * 86400), "s")))
        if len(later) + 1 < min_events:
            continue
        r = hypo_dist_km(d.latitude[i], d.longitude[i], d.depth[i], d.latitude.values[later],
                         d.longitude.values[later], d.depth.values[later])
        members = [i] + list(later[r <= km])
        if len(members) >= min_events:
            out.append(dict(first=i, members=members, ids=tuple(d.id[members])))
    return d, out


def resolve_overlaps(cands):
    """Drop candidates whose member set is a strict subset of another candidate's (documented rule)."""
    sets = [set(c["ids"]) for c in cands]
    return [c for k, c in enumerate(cands) if not any(sets[k] < s for s in sets)]


def migration_fit(t_days, r_km, n_perm=2000, seed=0):
    """Fit r = a + v t, r = v t and r = sqrt(4 pi D t) to events after the first (t > 0). Returns slopes,
    standard errors, residual RMS, envelope D, fitted D (m^2/s) and a permutation p-value of the r-t trend."""
    t, r = np.asarray(t_days, float), np.asarray(r_km, float)
    ok = t > 0
    t, r = t[ok], r[ok]
    res = dict(n_after_first=len(t))
    if len(t) < 2:
        return res
    ts, rm = t * 86400.0, r * 1000.0
    res["D_envelope_m2_s"] = float(np.max(rm ** 2 / (4 * np.pi * ts)))
    d_fit = float((rm ** 2 * ts).sum() / (4 * np.pi * (ts ** 2).sum()))
    res["D_fit_m2_s"] = d_fit
    res["rms_diffusive_km"] = float(np.sqrt(np.mean((np.sqrt(4 * np.pi * d_fit * ts) / 1000 - r) ** 2)))
    v0 = float((r * t).sum() / (t * t).sum())
    res["v_origin_km_day"] = v0
    res["rms_origin_km"] = float(np.sqrt(np.mean((v0 * t - r) ** 2)))
    if len(t) >= 3 and np.ptp(t) > 0:
        (v, a), cov = np.polyfit(t, r, 1, cov=True)
        res.update(v_free_km_day=float(v), v_free_se=float(np.sqrt(cov[0, 0])), intercept_km=float(a),
                   rms_free_km=float(np.sqrt(np.mean((v * t + a - r) ** 2))))
        rng = np.random.default_rng(seed)
        obs = abs(np.corrcoef(t, r)[0, 1])
        perm = [abs(np.corrcoef(t, rng.permutation(r))[0, 1]) for _ in range(n_perm)]
        res["trend_perm_p"] = float((np.sum(np.array(perm) >= obs) + 1) / (n_perm + 1))
    return res
