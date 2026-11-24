"""Checks on the real data files in data/ (read only). Skipped if the files are missing."""
import os

import numpy as np
import pandas as pd
import pytest

import seis_utils as su

need = pytest.mark.skipif(not (os.path.exists(su.COMBINED_CSV) and os.path.exists(su.USGS_CSV)),
                          reason="data files not present")


@need
def test_usgs_ids_unique_and_all_have_depth_type():
    us = pd.read_csv(su.USGS_CSV)
    dt = pd.read_csv(su.DEPTH_TYPES)
    assert us.id.is_unique
    assert set(us.id) == set(dt.id)
    assert us.mag.min() >= 3.0                 # the file name promises M3+


@need
def test_combined_keeps_every_usgs_event_once():
    us = pd.read_csv(su.USGS_CSV)
    co = pd.read_csv(su.COMBINED_CSV, dtype={"isc_id": str})
    assert co.id.is_unique
    assert set(us.id) <= set(co.id)
    assert co.isc_id.dropna().is_unique        # each ISC event used at most once


@need
def test_combined_matches_within_window():
    co = pd.read_csv(su.COMBINED_CSV, dtype={"isc_id": str})
    m = co[co.source == "USGS+ISC"]
    assert m.match_dt_s.notna().all() and (m.match_dt_s <= 20).all()
    assert (m.match_dist_km <= 100).all()
    assert co[co.source != "USGS+ISC"].match_dt_s.isna().all()


@need
def test_duplicate_flags_point_to_kept_rows():
    co = pd.read_csv(su.COMBINED_CSV, dtype={"isc_id": str})
    flagged = co[co.possible_duplicate_of.notna()]
    kept = set(co.loc[co.possible_duplicate_of.isna(), "id"])
    assert set(flagged.possible_duplicate_of) <= kept
    assert (flagged.source == "ISC").any()


@need
def test_load_catalog_counts_are_consistent():
    df, info = su.load_catalog("combined", include_fixed=False, start_year=1990, verbose=False)
    assert info["in_box_with_mag"] >= info["after_start_year"] >= info["after_duplicates"] >= info["after_depth_filter"]
    assert len(df) == info["after_depth_filter"]
    assert df.depth_type.isin(su.SOLVED_DEPTH).all()
    assert set(df.depth_class) <= {"intermediate", "shallow", "uncertain", "deep"}


@need
def test_isc_parser_reproduces_saved_events_for_one_year():
    import fetch_isc_bulletin as f
    path = os.path.join(f.RAW, "isc_isf_2013.txt")
    if not os.path.exists(path):
        pytest.skip("raw ISC file missing")
    rows = pd.DataFrame([r for r in map(f.summarise, f.parse(path)) if r])
    saved = pd.read_csv(f.OUT_EVENTS, dtype={"isc_id": str})
    saved = saved[saved.isc_id.isin(rows.isc_id)].set_index("isc_id").loc[rows.isc_id]
    assert len(saved) == len(rows)
    assert np.allclose(saved.depth.to_numpy(float), rows.depth.to_numpy(float), equal_nan=True)
    assert (saved.depth_type.to_numpy() == rows.depth_type.to_numpy()).all()
