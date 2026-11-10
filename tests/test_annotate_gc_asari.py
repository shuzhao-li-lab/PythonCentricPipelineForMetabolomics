"""Tests for ``pcpfm.annotate.gc_asari``. Run with `python -m unittest tests.test_annotate_gc_asari -v`."""

import os
import tempfile
import unittest
from unittest import mock

import pandas as pd

from pcpfm.annotate.gc_asari import run_asari_gc, write_asari_input
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

    def test_asari_reader_accepts_float_counts(self):
        from asari.tools.file_io import read_features_from_asari_table

        write_asari_input(gc_table(), ["S1", "S2"], self.path)
        with open(self.path) as f:
            n_samples, features = read_features_from_asari_table(f.read())
        self.assertEqual((n_samples, len(features)), (2, 3))
        self.assertEqual(features[0]["detection_counts"], 2)

    def test_missing_feature_value_fails(self):
        table = gc_table()
        table.loc[0, "snr"] = None
        with self.assertRaisesRegex(ValueError, "snr"):
            write_asari_input(table, ["S1", "S2"], self.path)

    def test_fractional_integer_column_fails(self):
        table = gc_table()
        table.loc[0, "detection_counts"] = 1.5
        with self.assertRaisesRegex(ValueError, "detection_counts"):
            write_asari_input(table, ["S1", "S2"], self.path)


PARAMS = {
    "ri_tolerance": 30,
    "score_cutoff_cosine": 0.4,
    "score_cutoff_entropy": 0.4,
    "corr_cutoff": 0.6,
    "denovo": False,
    "max_core_features": 20000,
}


def fake_gc_annotation(infile, outdir, KovatsIndex, database_file, **kwargs):
    """Stand-in for Asari: print like Asari and write a header-only Features table."""
    print("Done targeted annotation.")
    name = f"Features_{kwargs['project_name_handle']}.tsv"
    with open(os.path.join(outdir, name), "w") as f:
        f.write("feature\tmz\trtime\n")


class TestRunAsariGc(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.out = tmp.name

    def run_gc(self, params):
        return run_asari_gc(
            gc_table(), ["S1", "S2"], "k.tsv", "lib.msp", params, self.out
        )

    def test_passes_user_and_fixed_parameters(self):
        for denovo in (False, True):
            with mock.patch(
                "asari.annotate.annotate_gcms_full", side_effect=fake_gc_annotation
            ) as call:
                result = self.run_gc({**PARAMS, "denovo": denovo})
            kwargs = call.call_args.kwargs
            expected = {**PARAMS, "denovo": denovo, "max_ri_delta": 100}
            self.assertEqual(kwargs, kwargs | expected)
            self.assertFalse(kwargs["do_mirror_plot"])
            self.assertEqual(call.call_args.args[2:], ("k.tsv", "lib.msp"))
            self.assertTrue(os.path.exists(result["features_tsv"]))
            with open(result["log"]) as f:
                self.assertIn("Done targeted annotation", f.read())

    def test_failure_reports_log_tail(self):
        def broken(*args, **kwargs):
            print("Imported 0 compound library entries.")
            raise KeyError("ExactMass")

        with (
            mock.patch("asari.annotate.annotate_gcms_full", side_effect=broken),
            self.assertRaisesRegex(RuntimeError, "Imported 0 compound"),
        ):
            self.run_gc(PARAMS)

    def test_missing_features_table_fails(self):
        with (
            mock.patch("asari.annotate.annotate_gcms_full"),
            self.assertRaisesRegex(RuntimeError, "wrote no Features"),
        ):
            self.run_gc(PARAMS)
