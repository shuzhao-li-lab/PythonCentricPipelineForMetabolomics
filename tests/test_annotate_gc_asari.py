"""Tests for ``pcpfm.annotate.gc_asari``. Run with `python -m unittest tests.test_annotate_gc_asari -v`."""

import os
import tempfile
import unittest

import pandas as pd

from pcpfm.annotate.gc_asari import write_asari_input
from pcpfm.tables import FEATURE_COLUMNS


def gc_table(n=3, samples=("S1", "S2")):
    """A pcpfm feature table with Asari columns, samples, and a blank-filter column."""
    table = pd.DataFrame(
        {
            "id_number": [f"F{i}" for i in range(n)],
            "mz": [73.0 + i for i in range(n)],
            "rtime": [400.0 + 10 * i for i in range(n)],
            "rtime_left_base": [399.0 + 10 * i for i in range(n)],
            "rtime_right_base": [401.0 + 10 * i for i in range(n)],
            "parent_masstrack_id": list(range(n)),
            "peak_area": [1e6] * n,
            "cSelectivity": [1.0] * n,
            "goodness_fitting": [0.9] * n,
            "snr": [50] * n,
            "detection_counts": [float(len(samples))] * n,
            "blank_ratio": [10.0] * n,
        }
    )
    for k, s in enumerate(samples):
        table[s] = [1e5 * (k + 1) * (i + 1) for i in range(n)]
    return table


class TestAsariInput(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = os.path.join(tmp.name, "in.tsv")

    def test_samples_follow_the_11_feature_columns(self):
        table = gc_table()[["blank_ratio", "S2", *FEATURE_COLUMNS, "S1"]]
        write_asari_input(table, ["S1", "S2"], self.path)
        out = pd.read_csv(self.path, sep="\t")
        self.assertEqual(list(out.columns), FEATURE_COLUMNS + ["S1", "S2"])

    def test_missing_feature_column_fails(self):
        with self.assertRaisesRegex(ValueError, "snr"):
            write_asari_input(gc_table().drop(columns="snr"), ["S1", "S2"], self.path)
