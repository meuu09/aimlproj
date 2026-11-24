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
