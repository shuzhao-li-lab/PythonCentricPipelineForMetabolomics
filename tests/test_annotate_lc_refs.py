"""Tests for the reference loaders and searches of ``pcpfm annotate lc``.

Run with `python -m unittest tests.test_annotate_lc_refs -v`.
"""

import os
import tempfile
import unittest

import pandas as pd

from pcpfm.annotate.lc_evidence import EVIDENCE_COLUMNS
from pcpfm.annotate.lc_references import STANDARD_COLUMNS, load_standards
from pcpfm.annotate.lc_search import search_standards_ms1
from tests.test_annotate_lc import write_standards

EXAMPLES = os.path.join(os.path.dirname(__file__), "..", "examples", "auth_std_libraries")


class TestLoadStandards(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name

    def test_rt_in_seconds_and_first_precursor_column(self):
        csv = write_standards(os.path.join(self.d, "s.csv"), [("89675", "+", 1.5, 100.0)])
        std, warnings = load_standards(csv, "pos")
        self.assertEqual(list(std.columns), STANDARD_COLUMNS)
        row = std.iloc[0]
        self.assertEqual(row["rt_s"], 90.0)
        self.assertEqual(row["precursor_mz"], 100.0)
        # numeric compound names stay names
        self.assertEqual((row["reference_name"], row["reference_id"]), ("89675", "s.csv:2"))
        self.assertEqual(warnings, [])

    def test_other_polarity_dropped(self):
        csv = write_standards(
            os.path.join(self.d, "s.csv"), [("A", "+", 1.0, 100.0), ("B", "-", 1.0, 200.0)]
        )
        std, warnings = load_standards(csv, "neg")
        self.assertEqual(list(std["reference_name"]), ["B"])
        self.assertEqual(list(std["reference_id"]), ["s.csv:3"])
        self.assertIn("1 standards of the other polarity", warnings[0])
        self.assertEqual(len(load_standards(csv, None)[0]), 2)

    def test_rows_without_polarity_kept(self):
        csv = write_standards(
            os.path.join(self.d, "s.csv"), [("A", "", 1.0, 100.0), ("B", "-", 1.0, 200.0)]
        )
        std, warnings = load_standards(csv, "pos")
        self.assertEqual(list(std["reference_name"]), ["A"])
        self.assertIn("1 standards of the other polarity", warnings[0])

    def test_blank_lines_keep_line_numbers(self):
        csv = write_standards(os.path.join(self.d, "s.csv"), [("A", "+", 1.0, 100.0)])
        with open(csv, "a") as f:
            f.write("\nB,C1H1,+,1.0,200.0,50.0,201.0\n")
        std, _ = load_standards(csv, "pos")
        self.assertEqual(list(std["reference_id"]), ["s.csv:2", "s.csv:4"])

    def test_header_only(self):
        csv = write_standards(os.path.join(self.d, "s.csv"), [])
        with self.assertRaisesRegex(ValueError, "no standards"):
            load_standards(csv, "pos")

    def test_no_standard_for_polarity(self):
        csv = write_standards(os.path.join(self.d, "s.csv"), [("A", "+", 1.0, 100.0)])
        with self.assertRaisesRegex(ValueError, "no standards for neg"):
            load_standards(csv, "neg")

    def test_missing_column(self):
        csv = os.path.join(self.d, "s.csv")
        pd.DataFrame({"CompoundName": ["A"], "RT": [1.0]}).to_csv(csv, index=False)
        with self.assertRaisesRegex(ValueError, "Confirm Precursor"):
            load_standards(csv, "pos")

    def test_bad_value_names_line(self):
        csv = write_standards(
            os.path.join(self.d, "s.csv"), [("A", "+", 1.0, 100.0), ("B", "+", -1.0, 100.0)]
        )
        with self.assertRaisesRegex(ValueError, r"lines \[3\]"):
            load_standards(csv, "pos")

    def test_example_libraries(self):
        std, _ = load_standards(os.path.join(EXAMPLES, "auth_std_HILIC_pos_12_01_2024_JM_wMS2.csv"), "pos")
        self.assertEqual(len(std), 363)
        self.assertAlmostEqual(std["rt_s"].max(), 4.216290309 * 60)
        std, _ = load_standards(os.path.join(EXAMPLES, "auth_std_RP_neg_12_01_2024_JM_wMS2.csv"), "neg")
        self.assertEqual((len(std), set(std["polarity"])), (167, {"-"}))


def standards(rows):
    """Standards table from (name, m/z, RT seconds) tuples."""
    return pd.DataFrame(
        [
            {"reference_id": f"s.csv:{i + 2}", "reference_name": n, "neutral_formula": "C1",
             "precursor_mz": mz, "rt_s": rt, "polarity": "+", "reference_source": "s.csv"}
            for i, (n, mz, rt) in enumerate(rows)
        ],
        columns=STANDARD_COLUMNS,
    )


def features(rows):
    return pd.DataFrame(rows, columns=["id_number", "mz", "rtime"])


class TestSearchStandardsMs1(unittest.TestCase):
    def test_ppm_and_rt_boundaries(self):
        feats = features(
            [
                ("in", 100 * (1 + 4.9e-6), 69.9),
                ("ppm_out", 100 * (1 + 5.1e-6), 60.0),
                ("rt_out", 100.0, 70.1),
                # an RT in minutes must not pass for seconds
                ("minutes", 100.0, 1.0),
            ]
        )
        rows = search_standards_ms1(feats, standards([("A", 100.0, 60.0)]), None, 5, 10)
        self.assertEqual(list(rows["feature_id"]), ["in"])
        self.assertAlmostEqual(rows["mz_error_ppm"].iloc[0], 4.9, places=6)
        self.assertAlmostEqual(rows["rt_error_s"].iloc[0], 9.9)
        self.assertEqual(list(rows.columns[rows.columns.isin(EVIDENCE_COLUMNS)]), list(rows.columns))
        self.assertEqual((rows["level"].iloc[0], rows["evidence_type"].iloc[0]), ("1b", "standard_mz_rt"))

    def test_two_features_one_standard(self):
        feats = features([("F1", 100.0, 58.0), ("F2", 100.0001, 61.0)])
        rows = search_standards_ms1(feats, standards([("A", 100.0, 60.0)]), None, 5, 10)
        self.assertEqual(sorted(rows["feature_id"]), ["F1", "F2"])

    def test_closest_row_of_a_compound_kept(self):
        std = standards([("A", 100.0, 60.0), ("A", 100.0, 65.0), ("B", 100.0, 66.0)])
        rows = search_standards_ms1(features([("F1", 100.0, 64.0)]), std, None, 5, 10)
        self.assertEqual(sorted(rows["reference_name"]), ["A", "B"])
        a = rows[rows["reference_name"] == "A"].iloc[0]
        self.assertEqual((a["reference_id"], a["rt_error_s"]), ("s.csv:3", -1.0))

    def test_input_order_does_not_matter(self):
        feats = features([("F1", 100.0, 60.0), ("F2", 200.0, 30.0), ("F3", 100.0002, 55.0)])
        std = standards([("A", 100.0, 60.0), ("B", 200.0, 31.0)])
        key = ["feature_id", "reference_name"]
        a = search_standards_ms1(feats, std, None, 5, 10).sort_values(key).reset_index(drop=True)
        b = search_standards_ms1(feats.iloc[::-1], std.iloc[::-1], None, 5, 10)
        pd.testing.assert_frame_equal(a, b.sort_values(key).reset_index(drop=True))

    def test_empcpds_of_the_feature(self):
        membership = pd.DataFrame(
            {"feature_id": ["F1", "F1", "F2"], "empcpd_id": ["kp2", "_singleton_7", "kp3"]}
        )
        rows = search_standards_ms1(features([("F1", 100.0, 60.0)]), standards([("A", 100.0, 60.0)]), membership, 5, 10)
        self.assertEqual(rows["empcpd_id"].iloc[0], "_singleton_7;kp2")

    def test_no_match(self):
        rows = search_standards_ms1(features([("F1", 300.0, 60.0)]), standards([("A", 100.0, 60.0)]), None, 5, 10)
        self.assertEqual((len(rows), list(rows.columns)), (0, EVIDENCE_COLUMNS))


if __name__ == "__main__":
    unittest.main()
