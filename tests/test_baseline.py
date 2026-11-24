import os

import pandas as pd
import pytest

import seis_utils as su

BASE = os.path.join(su.HERE, "results_baseline")
needs_baseline = pytest.mark.skipif(not os.path.exists(os.path.join(BASE, "baseline_bvalues.csv")),
                                    reason="run reproduce_baseline.py first")


@pytest.fixture
def released():
    return pd.read_csv(os.path.join(BASE, "baseline_bvalues.csv")).set_index("group")


@needs_baseline
def test_released_code_full_catalog_matches_paper(released):
    assert released.loc["Full catalog", "N"] == 414
    assert released.loc["Full catalog", "Mc_released_code"] == pytest.approx(2.4)
    assert released.loc["Full catalog", "b_lsq"] == pytest.approx(0.91, abs=0.01)


@needs_baseline
def test_released_code_cluster_values_are_stable(released):
    assert released.loc["Cluster 1", "N"] == 254
    assert released.loc["Cluster 1", "Mc_released_code"] == pytest.approx(2.4)
    assert released.loc["Cluster 2", "N"] == 61


@needs_baseline
def test_cluster1_b_at_published_mc_agrees_within_sampling_error():
    m = pd.read_csv(os.path.join(BASE, "baseline_magnitudes.csv"))
    c1 = m[m.group == "Cluster 1"].mag.to_numpy()
    b, sigma, n = su.b_mle(c1, 2.85, dm=0.01)
    assert abs(b - 0.96) <= sigma and n >= 50


@needs_baseline
@pytest.mark.xfail(strict=True, reason="Published Cluster 1 Mc 2.85 is not produced by the released code (0.5 bins)")
def test_published_cluster1_mc(released):
    assert released.loc["Cluster 1", "Mc_released_code"] == pytest.approx(2.85)


@needs_baseline
@pytest.mark.xfail(strict=True, reason="Published Cluster 1 has 248 events; the released polygons select 254")
def test_published_cluster1_count(released):
    assert released.loc["Cluster 1", "N"] == 248
