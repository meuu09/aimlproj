"""Unit tests for seis_utils. Inputs are tiny hand-made tables that exercise the code;
they are test fixtures, not earthquake data, and are never written anywhere."""
import numpy as np
import pandas as pd
import pytest

import seis_utils as su


def test_max_curv_puts_magnitudes_in_their_own_bin():
    assert su.max_curv([4.4] * 5 + [4.5] * 2 + [4.2]) == 4.4
    assert su.max_curv([3.0, 3.1, 3.1, 3.2]) == 3.1


def test_max_curv_empty_is_nan():
    assert np.isnan(su.max_curv([]))


def _eqcorrscan_verbatim(magnitudes, bin_size=0.5):
    min_bin, max_bin = int(min(magnitudes)), int(max(magnitudes) + 1)
    bins = np.arange(min_bin, max_bin + bin_size, bin_size)
    df, bins = np.histogram(magnitudes, bins)
    grad = (df[1:] - df[0:-1]) / bin_size
    curvature = (grad[1:] - grad[0:-1]) / bin_size
    return float(bins[np.argmax(np.abs(curvature))] + bin_size)


@pytest.mark.parametrize("mags", [
    [3.1, 3.4, 3.6, 3.9, 4.0, 4.0, 4.1, 4.2, 4.4, 4.6, 5.0, 5.3],
    [1.2, 2.0, 2.2, 2.5, 2.5, 2.6, 2.7, 3.0, 3.3, 4.1],
])
def test_max_curv_eqcorrscan_matches_original(mags):
    assert su.max_curv_eqcorrscan(mags) == _eqcorrscan_verbatim(mags)


def test_b_mle_formula_and_error():
    mags = [4.4] * 10
    b, err, n = su.b_mle(mags, 4.4)
    assert n == 10
    assert b == pytest.approx(np.log10(np.e) / 0.05)
    assert err == pytest.approx(b / np.sqrt(10))


def test_b_mle_includes_events_at_float_threshold():
    mc = 4.0 + 0.4
    assert su.b_mle([4.4] * 6, mc)[2] == 6


def test_b_mle_too_few_events():
    b, err, n = su.b_mle([4.5, 4.6, 4.7], 4.4)
    assert np.isnan(b) and np.isnan(err) and n == 3


def test_b_lsq_exact_gutenberg_richter():
    mags = [3.0] * 900 + [4.0] * 90 + [5.0] * 10
    b, steps, cum, fn = su.b_lsq(mags, 3.0)
    assert b == pytest.approx(1.0)
    assert list(cum) == [1000, 100, 10]


def test_b_lsq_too_few_steps_returns_nan():
    b, steps, cum, fn = su.b_lsq([4.5, 4.5, 4.6], 4.4)
    assert np.isnan(b) and fn is None


@pytest.mark.parametrize("depth,err,expected", [
    (50, 10, "intermediate"),
    (49, 10, "uncertain"),
    (30, 9.9, "shallow"),
    (30, 10, "uncertain"),
    (60, np.nan, "uncertain"),
    (145, 5, "intermediate"),
    (150, 5, "uncertain"),
    (151, 5, "uncertain"),
    (156, 5, "deep"),
    (210, np.nan, "uncertain"),
    (10, 0, "shallow"),
])
def test_depth_class(depth, err, expected):
    assert su.depth_class([depth], [err])[0] == expected


def test_interevent_days_sorted_and_cv():
    t = pd.Series(pd.to_datetime(["2000-01-03", "2000-01-01", "2000-01-02"]))
    _, gaps = su.interevent_days(t)
    assert list(gaps) == [1.0, 1.0]
    assert su.cv(gaps) == 0.0
    assert np.isnan(su.cv([5.0]))


def test_haversine_one_degree_latitude():
    assert su.haversine_km(33.0, 76.0, 34.0, 76.0) == pytest.approx(111.19, abs=0.05)


def _cat(rows):
    return pd.DataFrame(rows, columns=["time", "latitude", "longitude"])


def test_match_events_one_to_one_closest_first():
    a = _cat([("2010-01-01T00:00:00Z", 33.0, 76.0),
              ("2010-01-01T00:00:05Z", 33.1, 76.0)])
    b = _cat([("2010-01-01T00:00:01", 33.0, 76.0),
              ("2010-01-01T00:00:05", 33.1, 76.0),
              ("2010-01-01T00:00:30", 33.0, 76.0),
              ("2010-01-01T00:00:00", 35.0, 76.0)])
    m = su.match_events(a, b, 20, 100)
    assert sorted(zip(m.ia, m.ib)) == [(0, 0), (1, 1)]


def test_match_events_each_event_used_once():
    a = _cat([("2010-01-01T00:00:00Z", 33.0, 76.0), ("2010-01-01T00:00:01Z", 33.0, 76.0)])
    b = _cat([("2010-01-01T00:00:00", 33.0, 76.0)])
    m = su.match_events(a, b, 20, 100)
    assert len(m) == 1 and m.ia[0] == 0


def _dup_cat():
    return pd.DataFrame({
        "id": ["us1", "isc2", "us3", "us4", "isc5"],
        "time": ["2010-01-01T00:00:00Z", "2010-01-01T00:00:30", "2011-01-01T00:00:00Z",
                 "2011-01-01T00:00:05Z", "2012-01-01T00:00:10"],
        "latitude": [33.0, 33.05, 33.0, 33.0, 33.0],
        "longitude": [76.0, 76.0, 76.0, 76.1, 76.0],
        "source": ["USGS+ISC", "ISC", "USGS", "USGS", "ISC"],
        "depth_type": ["from location"] * 5,
        "depthError": [5.0] * 5,
        "mag": [4.5, 4.3, 4.0, 4.1, 3.0],
    })


def test_flag_duplicates_rules():
    c = _dup_cat()
    dup, pairs = su.flag_duplicates(c)
    assert dup.loc[1] == "us1"
    assert pd.isna(dup.loc[0])
    assert pd.isna(dup.loc[2]) and pd.isna(dup.loc[3])
    assert pd.isna(dup.loc[4])
    assert set(pairs.flagged) == {True, False}


def test_flag_duplicates_ignores_far_pairs():
    c = _dup_cat()
    c.loc[1, "latitude"] = 33.6
    dup, _ = su.flag_duplicates(c)
    assert dup.isna().all()


def test_parse_args_rejects_drop_duplicates_for_usgs():
    with pytest.raises(SystemExit):
        su.parse_args(argv=["--drop-duplicates"])


def test_parse_args_rejects_unknown_catalog():
    with pytest.raises(SystemExit):
        su.parse_args(argv=["--catalog", "isc"])


def test_parse_args_combined_nodup_folder():
    a = su.parse_args(argv=["--catalog", "combined", "--drop-duplicates", "--include-fixed"])
    assert su.results_base(a).endswith("combined_usgs_isc_nodup")
    assert su.case_name(a) == "including_fixed"


def test_shi_bolt_sigma_formula():
    m = np.array([4.4, 4.5, 4.6, 4.8, 5.0, 5.3])
    b = su.b_mle(m, 4.4)[0]
    expected = 2.30 * b ** 2 * np.sqrt(((m - m.mean()) ** 2).sum() / (len(m) * (len(m) - 1)))
    assert su.b_shi_bolt_sigma(m, 4.4) == pytest.approx(expected)
    assert np.isnan(su.b_shi_bolt_sigma([4.5, 4.6], 4.4))


def test_bootstrap_is_repeatable_and_brackets_mle():
    m = [4.4, 4.4, 4.5, 4.5, 4.6, 4.7, 4.8, 5.0, 5.2, 5.6]
    a, b = su.b_bootstrap(m, 4.4, seed=1), su.b_bootstrap(m, 4.4, seed=1)
    assert a == b
    med, lo, hi = a
    assert lo < su.b_mle(m, 4.4)[0] < hi
    assert all(np.isnan(su.b_bootstrap([4.5], 4.4)))


def test_mc_sliding_windows():
    t = pd.date_range("2000-01-01", periods=30, freq="D")
    mags = np.r_[[3.5] * 15, [4.5] * 15]
    r = su.mc_sliding(t, mags, window=10, step=10, n_boot=20)
    assert len(r) == 3
    assert list(r.mc_maxc) == [3.5, 3.5, 4.5]
    assert (r.mc_maxc_p16 <= r.mc_maxc_p84).all()


def test_cluster_outline_contains_kishtwar():
    from matplotlib.path import Path
    poly = Path(su.CLUSTER_POLY)
    assert poly.contains_point((75.8, 33.3))
    assert not poly.contains_point((77.5, 32.0))


def test_mc_sliding_narrow_window_gives_nan_for_authors_method():
    r = su.mc_sliding(pd.date_range("2000-01-01", periods=10, freq="D"), [3.5] * 10, window=10, step=10, n_boot=5)
    assert np.isnan(r.mc_authors.iloc[0]) and r.mc_maxc.iloc[0] == 3.5


def test_duplicate_chain_never_joins_two_usgs_events():
    c = pd.DataFrame({
        "id": ["usA", "iscX", "usB"],
        "time": ["2010-01-01T00:00:00Z", "2010-01-01T00:00:40", "2010-01-01T00:01:20Z"],
        "latitude": [33.0, 33.0, 33.0], "longitude": [76.0, 76.0, 76.0],
        "source": ["USGS", "ISC", "USGS"], "depth_type": ["from location"] * 3,
        "depthError": [5.0] * 3, "mag": [4.5, 4.4, 4.3]})
    dup, pairs = su.flag_duplicates(c)
    assert pd.isna(dup.loc[0]) and pd.isna(dup.loc[2])
    assert dup.loc[1] == "usA"
    assert pairs.flagged.sum() == 1


def test_duplicate_isc_only_group_keeps_best_row():
    c = pd.DataFrame({
        "id": ["isc1", "isc2"], "time": ["2010-01-01T00:00:00", "2010-01-01T00:00:20"],
        "latitude": [33.0, 33.0], "longitude": [76.0, 76.0], "source": ["ISC", "ISC"],
        "depth_type": ["unknown", "from location"], "depthError": [np.nan, 5.0], "mag": [4.0, 3.9]})
    dup, _ = su.flag_duplicates(c)
    assert dup.loc[0] == "isc2" and pd.isna(dup.loc[1])


def test_pair_flag_is_pair_level():
    c = _dup_cat()
    _, pairs = su.flag_duplicates(c)
    p = pairs.set_index(["id_1", "id_2"])
    assert p.loc[("us1", "isc2"), "flagged"]
    assert not p.loc[("us3", "us4"), "eligible"] and not p.loc[("us3", "us4"), "flagged"]


def test_match_events_competing_candidates_resolved_by_score():
    a = _cat([("2010-01-01T00:00:00Z", 33.0, 76.0), ("2010-01-01T00:00:10Z", 33.0, 76.0)])
    b = _cat([("2010-01-01T00:00:09", 33.0, 76.0)])
    m = su.match_events(a, b, 20, 100)
    assert len(m) == 1 and m.ia[0] == 1


def test_norm_magtype_keeps_scales_apart():
    out = list(su.norm_magtype(["MW", "mww", "mwc", "mb", "ML", "mL", "MD", "Ms", None, "Mx"]))
    assert out == ["Mw", "Mw", "Mw", "mb", "ML", "ML", "MD", "Ms", "none", "Mx"]


def test_load_catalog_respects_dates_and_box():
    df, _ = su.load_catalog("combined", include_fixed=True, verbose=False)
    assert df.time.min() >= su.CATALOG_START and df.time.max() < su.CATALOG_END
    assert df.latitude.between(*su.BOX["lat"]).all() and df.longitude.between(*su.BOX["lon"]).all()


def _synthetic_gr(n=4000, b=1.0, mc=3.0, seed=1):
    """Synthetic (clearly labelled) Gutenberg-Richter magnitudes, complete above mc, 0.1 binned,
    with a fading tail below mc."""
    rng = np.random.default_rng(seed)
    m = mc - 0.05 + rng.exponential(np.log10(np.e) / b, n)
    tail = mc - 0.05 - rng.exponential(0.3, n // 3)
    keep = rng.random(tail.size) < np.exp(-(mc - tail) * 2)
    return np.round(np.r_[m, tail[keep]], 1)


def test_gft_and_mbs_find_synthetic_completeness():
    m = _synthetic_gr()
    assert su.mc_gft(m, level=90) == pytest.approx(3.0, abs=0.2)
    assert su.mc_mbs(m) == pytest.approx(3.0, abs=0.3)
    assert np.isnan(su.mc_gft([3.0] * 5))


def test_b_mle_recovers_synthetic_b():
    m = _synthetic_gr(n=20000, b=1.0)
    b, s, n = su.b_mle(m, 3.0)
    assert b == pytest.approx(1.0, abs=3 * s)


def test_poisson_ci_known_values():
    lo, hi = su.poisson_ci(0, 1.0)
    assert lo == 0 and hi == pytest.approx(3.689, abs=1e-3)
    lo, hi = su.poisson_ci(10, 2.0)
    assert lo == pytest.approx(4.795 / 2, abs=1e-3) and hi == pytest.approx(18.390 / 2, abs=1e-3)


def test_cv_poisson_test_regular_series_is_low_and_clustered_is_high():
    t = pd.date_range("2000-01-01", periods=30, freq="30D")
    cv, lo, hi, p = su.cv_poisson_test(t, "2000-01-01", t[-1] + pd.Timedelta(days=30), n_sim=500)
    assert cv < lo and p < 0.05
    burst = list(pd.date_range("2000-01-01", periods=20, freq="h")) + list(pd.date_range("2005-01-01", periods=5, freq="365D"))
    cv, lo, hi, p = su.cv_poisson_test(burst, "2000-01-01", "2010-01-01", n_sim=500)
    assert cv > hi and p < 0.05
    assert np.isnan(su.cv_poisson_test(t[:3], "2000", "2001")[0])


def test_b_reliability_labels():
    assert su.b_reliability(50).startswith("reliable")
    assert su.b_reliability(20).startswith("rough")
    assert su.b_reliability(19).startswith("not reported")


def test_dispersion_index_constant_counts_is_zero():
    t = pd.to_datetime([f"{y}-06-01" for y in range(2000, 2010)])
    d, p = su.dispersion_index(t, "2000-01-01", "2009-12-31")
    assert d == 0 and p > 0.99


def test_hypo_dist_pythagoras():
    assert su.hypo_dist_km(33.0, 76.0, 50.0, 33.0, 76.0, 80.0) == pytest.approx(30.0)
    h = su.haversine_km(33.0, 76.0, 33.1, 76.0)
    assert su.hypo_dist_km(33.0, 76.0, 50.0, 33.1, 76.0, 60.0) == pytest.approx(np.hypot(h, 10.0))


def _seq_cat():
    return pd.DataFrame({"id": list("abcde"),
                         "time": pd.to_datetime(["2000-01-01", "2000-01-05", "2000-01-10", "2000-01-25", "2001-01-01"]),
                         "latitude": [33.0, 33.05, 33.1, 33.1, 33.0], "longitude": [76.0] * 5,
                         "depth": [60.0, 60.0, 60.0, 60.0, 60.0]})


def test_candidate_sequences_keeps_all_and_subset_rule():
    _, c = su.candidate_sequences(_seq_cat(), 20, 20, 3)
    assert [x["ids"] for x in c] == [("a", "b", "c"), ("b", "c", "d")]
    _, c2 = su.candidate_sequences(_seq_cat(), 30, 20, 2)
    kept = su.resolve_overlaps(c2)
    assert ("b", "c", "d") not in [x["ids"] for x in kept] and ("a", "b", "c", "d") in [x["ids"] for x in kept]


def test_candidate_sequences_window_boundaries_inclusive():
    c = _seq_cat()
    _, s = su.candidate_sequences(c.iloc[[0, 1]].reset_index(drop=True), 4, 20, 2)
    assert s and s[0]["ids"] == ("a", "b")
    _, s = su.candidate_sequences(c.iloc[[0, 1]].reset_index(drop=True), 3.99, 20, 2)
    assert not s


def test_migration_fit_recovers_linear_and_diffusive():
    t = np.array([0, 1, 2, 3, 4, 5.0])
    fit = su.migration_fit(t, 0.5 * t, n_perm=200)
    assert fit["v_origin_km_day"] == pytest.approx(0.5) and fit["v_free_km_day"] == pytest.approx(0.5)
    assert fit["trend_perm_p"] < 0.05
    D = 0.5
    r = np.sqrt(4 * np.pi * D * t * 86400) / 1000
    fit = su.migration_fit(t, r, n_perm=10)
    assert fit["D_fit_m2_s"] == pytest.approx(D) and fit["D_envelope_m2_s"] == pytest.approx(D)
    assert su.migration_fit([0, 1], [0, 1]) == {"n_after_first": 1}


def test_edge_cases_return_nan_not_crash():
    assert np.isnan(su.b_mle([], 4.4)[0]) and su.b_mle([], 4.4)[2] == 0
    assert su.b_mle([4.5] * 10, 4.5)[0] == pytest.approx(np.log10(np.e) / 0.05)
    assert np.isnan(su.b_lsq([4.5] * 10, 4.4)[0])
    assert np.isnan(su.b_bootstrap([4.5, 4.6], 4.4)[0])
    assert np.isnan(su.cv([])) and np.isnan(su.cv([0.0, 0.0]))
    t, g = su.interevent_days(pd.Series([], dtype="datetime64[ns]"))
    assert len(g) == 0
    assert su.depth_class([np.nan], [5.0])[0] == "uncertain"
    assert np.isnan(su.max_curv([np.nan]))
    assert np.isnan(su.mc_mbs([4.0] * 10)) and np.isnan(su.mc_gft([]))
    assert np.isnan(su._eqcorrscan_or_nan([4.0, 4.1]))
    d, c = su.candidate_sequences(_seq_cat().iloc[:0], 20, 20, 3)
    assert c == []
    fit = su.migration_fit([0, 1, 1, 1], [0, 2, 2, 2], n_perm=10)
    assert "v_origin_km_day" in fit and np.isfinite(fit["v_origin_km_day"])


def test_load_catalog_unknown_catalog_raises():
    with pytest.raises(ValueError):
        su.load_catalog("isc", verbose=False)


GENUINE_ISF = """Event     8443952 Kashmir-India border region
   Date       Time        Err   RMS Latitude Longitude  Smaj  Smin  Az Depth   Err Ndef Nsta Gap  mdist  Mdist Qual   Author      OrigID
2005/01/04 21:13:51.50   3.23 1.000  33.9890   76.4380  14.5 25.60  -1  19.7  14.7        15 306                      NDI        7392851

Magnitude  Err Nsta Author      OrigID
MD     3.8          NDI        7392851
ML     3.8          NDI        7392851

Event     8500699 Western Xizang-India border region
   Date       Time        Err   RMS Latitude Longitude  Smaj  Smin  Az Depth   Err Ndef Nsta Gap  mdist  Mdist Qual   Author      OrigID
2005/02/02 14:48:42.00   3.98 1.300  31.7050   78.5430 10.10 11.30  -1   5.0f             15 215                      NDI        7433615
"""


def test_isf_parser_on_genuine_records(tmp_path):
    """Genuine ISC Bulletin lines (data/isc_raw/isc_isf_2005.txt), values checked by hand."""
    from fetch import fetch_isc_bulletin as fib
    p = tmp_path / "isf.txt"
    p.write_text(GENUINE_ISF)
    ev = list(fib.parse(str(p)))
    assert [e["id"] for e in ev] == ["8443952", "8500699"]
    o = ev[0]["origins"][0]
    assert o["time"] == "2005-01-04T21:13:51.50" and o["time_err"] == 3.23 and o["rms"] == 1.0
    assert (o["lat"], o["lon"], o["smaj"], o["smin"]) == (33.989, 76.438, 14.5, 25.6)
    assert (o["depth"], o["depth_fixed"], o["depth_err"], o["ndef"], o["nsta"], o["gap"]) == (19.7, False, 14.7, None, 15, 306)
    assert (o["author"], o["origid"]) == ("NDI", "7392851")
    assert ev[0]["mags"] == [dict(type="MD", mag=3.8, author="NDI"), dict(type="ML", mag=3.8, author="NDI")]
    f = ev[1]["origins"][0]
    assert f["depth"] == 5.0 and f["depth_fixed"] and f["depth_err"] is None
    row = fib.summarise(ev[1])
    assert row["depth_type"] == "operator assigned" and row["depthError"] is None and row["origin_id"] == "7433615"
    assert fib.summarise(ev[0])["depth_type"] == "from location"


def test_dbscan_haversine_lat_lon_order():
    from robustness import cluster_robustness as cr
    lat = np.r_[np.full(5, 33.0) + np.arange(5) * 0.01, np.full(5, 33.0)]
    lon = np.r_[np.full(5, 76.0), np.full(5, 76.0 + 60 / 93.3) + np.arange(5) * 0.01]
    lab = cr.dbscan(lat, lon, 30, 4)
    assert len(set(lab[:5])) == 1 and len(set(lab[5:])) == 1 and lab[0] != lab[5] and (lab >= 0).all()
    assert (cr.dbscan(lat[:3], lon[:3], 30, 4) == -1).all()


def test_pca_azimuth_of_a_known_line():
    from analysis import clustering_checks as cc
    t = np.linspace(-50, 50, 21)
    az = np.radians(130.0)
    x, y = t * np.sin(az), t * np.cos(az)
    lon = cc.LON0 + x / (111.32 * np.cos(np.radians(cc.LAT0)))
    lat = cc.LAT0 + y / 110.57
    a, share = cc.pca_azimuth(pd.DataFrame({"longitude": lon, "latitude": lat}))
    assert a == pytest.approx(130.0, abs=0.5) and share == pytest.approx(1.0)


def test_monte_carlo_is_reproducible_and_skips_missing_errors():
    from robustness import cluster_robustness as cr
    df = pd.DataFrame({"id": [f"e{i}" for i in range(8)], "time": pd.date_range("2000", periods=8, freq="D"),
                       "latitude": 33.0 + np.arange(8) * 0.02, "longitude": [75.9] * 8,
                       "depth": [60.0] * 7 + [60.0], "depthError": [5.0] * 7 + [np.nan],
                       "horizontalError": [5.0] * 8, "mag": [4.5] * 8, "magType": ["mb"] * 8})
    df["depth_class"] = su.depth_class(df.depth, df.depthError)
    e1, r1 = cr.monte_carlo("t", df, 8.0, n_runs=20, seed=3)
    e2, r2 = cr.monte_carlo("t", df, 8.0, n_runs=20, seed=3)
    pd.testing.assert_frame_equal(r1, r2)
    assert np.isnan(e1.p_intermediate.iloc[-1]) and e1.depth_error_missing.iloc[-1]
    assert (e1.p_intermediate.iloc[:-1] == 1.0).all()
