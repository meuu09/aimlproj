import unittest

from catalog_qc import candidate_review_reason, depth_class, find_matches, haversine_km, parse_utc
from fetch_isc_bulletin import summarise
from fetch_usgs_catalog import query_url


class CatalogQcTests(unittest.TestCase):
    def test_parse_utc_treats_naive_isc_time_as_utc_and_truncates_ms(self):
        self.assertEqual(
            parse_utc("2020-01-02T03:04:05.123987"),
            parse_utc("2020-01-02T03:04:05.123Z"),
        )

    def test_haversine_same_point_is_zero(self):
        self.assertAlmostEqual(haversine_km(33.0, 76.0, 33.0, 76.0), 0.0)

    def test_greedy_matching_is_one_to_one_and_reports_ambiguity(self):
        usgs = [
            {"id": "usgs-a", "time": "2020-01-01T00:00:00.000Z", "latitude": 33, "longitude": 76},
            {"id": "usgs-b", "time": "2020-01-01T00:00:00.500Z", "latitude": 33, "longitude": 76.01},
        ]
        isc = [
            {"isc_id": "isc-a", "time": "2020-01-01T00:00:00.000", "latitude": 33, "longitude": 76}
        ]
        matches, candidates = find_matches(usgs, isc)
        self.assertEqual(candidates, 2)
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0]["usgs_index"], 0)
        self.assertEqual(matches[0]["isc_candidate_count"], 2)
        self.assertIn("multiple_usgs_candidates_for_isc", candidate_review_reason(matches[0]))

    def test_spatial_residual_is_review_only(self):
        usgs = [{"id": "u", "time": "2020-01-01T00:00:00Z", "latitude": 33, "longitude": 76}]
        isc = [{"isc_id": "i", "time": "2020-01-01T00:00:00Z", "latitude": 33, "longitude": 76.3}]
        matches, _ = find_matches(usgs, isc)
        self.assertEqual(len(matches), 1)  # within the retained 100 km match limit
        self.assertGreater(matches[0]["horizontal_distance_km"], 20)
        self.assertIn("horizontal_residual_over_review_distance", candidate_review_reason(matches[0]))

    def test_depth_class_keeps_large_and_missing_errors_uncertain(self):
        self.assertEqual(depth_class(55.0, 1000.0), "uncertain")
        self.assertEqual(depth_class(55.0, None), "uncertain")
        self.assertEqual(depth_class(52.0, 5.0), "intermediate")
        self.assertEqual(depth_class(20.0, 5.0), "shallow")

    def test_isc_phase_depth_flag_and_1000_km_error_are_preserved(self):
        event = {
            "id": "test-only",
            "region": "test fixture",
            "origins": [
                {
                    "time": "2020-01-01T00:00:00.00",
                    "lat": 33.0,
                    "lon": 76.0,
                    "depth": 52.0,
                    "depth_flag": "d",
                    "depth_fixed": False,
                    "depth_err": 1000.0,
                    "author": "NDI",
                    "origid": "phase-origin",
                    "prime": True,
                }
            ],
            "mags": [
                {"type": "mL", "mag": 3.4, "error": 0.2, "nsta": 4, "author": "NDI", "origid": "m1"}
            ],
        }
        row = summarise(event)
        self.assertEqual(row["depth_type"], "constrained by depth phases")
        self.assertEqual(row["depth_flag"], "d")
        self.assertEqual(row["depthError"], 1000.0)
        self.assertEqual(row["depth_error_confidence"], "90% (ISF)")
        self.assertTrue(row["depth_error_gt_100_km"])
        self.assertEqual(row["magType"], "mL")
        self.assertEqual(row["magError"], 0.2)

    def test_usgs_query_url_records_explicit_bounds_and_magnitude_floor(self):
        url = query_url({"format": "csv", "minmagnitude": "3.0", "starttime": "1970-03-05T00:00:00"})
        self.assertIn("format=csv", url)
        self.assertIn("minmagnitude=3.0", url)
        self.assertIn("starttime=1970-03-05T00%3A00%3A00", url)


if __name__ == "__main__":
    unittest.main()
