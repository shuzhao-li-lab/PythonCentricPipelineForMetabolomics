"""Explicit sample-column schemas survive annotation additions and QC filtering."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from pcpfm import cli, process, qc, tables
from pcpfm.project import ProjectRecord


class TestTables(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.record = ProjectRecord.load_or_create(str(self.root))
        self.ft = pd.DataFrame({
            "id_number": ["001", "002"],
            "mz": [100.0, 200.0],
            "rtime": [12.5, 24.0],
            "rtime_left_base": [12.0, 23.5],
            "rtime_right_base": [13.0, 24.5],
            "parent_masstrack_id": [1, 2],
            "peak_area": [1000.0, 2000.0],
            "cSelectivity": [0.9, 0.8],
            "goodness_fitting": [0.95, 0.92],
            "snr": [10.0, 20.0],
            "detection_counts": [2, 2],
            "new_annotation": ["A", "B"],
            "new_score": [0.1, 0.2],
            "S1": [10.0, 20.0],
            "S2": [30.0, 40.0],
        })
        self.path = self.root / "table.tsv"
        self.ft.to_csv(self.path, sep="\t", index=False)

    def register(self, **extra):
        self.record.add_artifact("table", str(self.path), "test table", "process", **extra)

    def test_schema_roundtrip_preserves_unknown_annotations_and_sample_order(self):
        self.register(sample_columns=["S2", "S1"])
        self.record.save()
        record = ProjectRecord.load(str(self.root))
        ft, cols = tables.load_table(record, "table")
        self.assertEqual(cols, ["S2", "S1"])
        pd.testing.assert_frame_equal(ft, self.ft)

    def test_legacy_uses_inputs_not_exclusion(self):
        self.register()
        self.record.data["inputs"] = [{"sample": "S1"}, {"sample": "S2"}, {"sample": "dropped"}]
        self.assertEqual(tables.load_table(self.record, "table")[1], ["S1", "S2"])
        self.record.data["inputs"] = []
        with self.assertRaisesRegex(ValueError, "no sample_columns"):
            tables.load_table(self.record, "table")

    def test_invalid_schema_does_not_fall_back(self):
        for cols in (["missing"], ["S1", "S1"], []):
            with self.subTest(columns=cols), self.assertRaises(ValueError):
                self.register(sample_columns=cols)
                tables.load_table(self.record, "table")

        self.register(sample_columns="S1")
        with self.assertRaisesRegex(TypeError, "sample_columns must be a list"):
            tables.load_table(self.record, "table")

    def test_duplicate_sample_columns_are_rejected(self):
        self.register(sample_columns=["S1", "S1"])

        with self.assertRaisesRegex(ValueError, "unique sample columns"):
            tables.load_table(self.record, "table")

    def test_table_without_sample_columns_is_rejected(self):
        path = self.root / "no_samples.tsv"
        pd.DataFrame({"id_number": ["001"], "mz": [100.0], "rtime": [12.5]}).to_csv(
            path, sep="\t", index=False
        )
        self.record.add_artifact("no_samples", str(path), "test table", "process")

        with self.assertRaisesRegex(ValueError, "no sample_columns"):
            tables.load_table(self.record, "no_samples")

    def test_table_without_feature_metadata_still_loads(self):
        path = self.root / "samples_only.tsv"
        samples_only = pd.DataFrame({"S1": [10.0, 20.0], "S2": [30.0, 40.0]})
        samples_only.to_csv(path, sep="\t", index=False)
        self.record.add_artifact(
            "samples_only", path, "sample-only table", "process",
            sample_columns=["S1", "S2"],
        )

        loaded, samples = tables.load_table(self.record, "samples_only")

        self.assertEqual(samples, ["S1", "S2"])
        pd.testing.assert_frame_equal(loaded, samples_only)

    def test_asari_uses_only_known_stems_and_allows_dropped_samples(self):
        ft, cols = tables.load_asari_table(self.path, {"S1": "renamed", "S2": "second", "absent": "dropped"})
        self.assertEqual(cols, ["renamed", "second"])
        self.assertEqual(ft["new_annotation"].tolist(), ["A", "B"])
        self.assertEqual(ft["new_score"].tolist(), [0.1, 0.2])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            tables.load_asari_table(self.path, {"S1": "new_annotation"})

    def test_expected_asari_table_loads_and_keeps_annotations(self):
        loaded, samples = tables.load_asari_table(
            self.path, {"S1": "sample_1", "S2": "sample_2"}
        )

        self.assertEqual(samples, ["sample_1", "sample_2"])
        self.assertTrue(set(tables.FEATURE_COLUMNS).issubset(loaded.columns))
        self.assertEqual(loaded["id_number"].tolist(), ["001", "002"])
        self.assertEqual(loaded["mz"].tolist(), [100.0, 200.0])
        self.assertEqual(loaded["new_annotation"].tolist(), ["A", "B"])
        self.assertEqual(loaded["new_score"].tolist(), [0.1, 0.2])
        self.assertEqual(loaded["sample_1"].tolist(), [10.0, 20.0])
        self.assertEqual(loaded["sample_2"].tolist(), [30.0, 40.0])

    def test_gibberish_table_is_rejected(self):
        path = self.root / "gibberish.tsv"
        pd.DataFrame({"nonsense": ["not a feature table"], "garbage": ["???"]}).to_csv(
            path, sep="\t", index=False
        )

        with self.assertRaisesRegex(ValueError, "non-empty list"):
            tables.load_asari_table(path, {"S1": "sample_1"})

    def test_all_zero_sample_is_loaded_without_being_dropped(self):
        zero_table = self.ft.copy()
        zero_table["S1"] = 0.0
        zero_table["S2"] = 0.0
        path = self.root / "all_zero.tsv"
        zero_table.to_csv(path, sep="\t", index=False)

        loaded, samples = tables.load_asari_table(
            path, {"S1": "sample_1", "S2": "sample_2"}
        )

        self.assertEqual(samples, ["sample_1", "sample_2"])
        self.assertTrue((loaded[["sample_1", "sample_2"]] == 0.0).all().all())

    def test_load_table_preserves_all_zero_sample(self):
        zero_table = self.ft.copy()
        zero_table["S1"] = 0.0
        zero_table["S2"] = 0.0
        zero_table.to_csv(self.path, sep="\t", index=False)
        self.register(sample_columns=["S1", "S2"])

        loaded, samples = tables.load_table(self.record, "table")

        self.assertEqual(samples, ["S1", "S2"])
        self.assertTrue((loaded[samples] == 0.0).all().all())

    def test_metadata_alignment_and_missing_metadata(self):
        info = tables.sample_info(self.record, ["S2", "S1"])
        self.assertEqual(info.index.tolist(), ["S2", "S1"])
        self.assertEqual(info["injection_order"].tolist(), [0, 1])
        path = self.root / "metadata.tsv"
        pd.DataFrame({"sample": ["S1", "S2"], "sample_type": ["sample", "pooled"]}).to_csv(path, sep="\t", index=False)
        self.record.add_artifact("metadata", str(path), "metadata", "process")
        info = tables.sample_info(self.record, ["S2", "S1"])
        self.assertEqual(info["injection_order"].tolist(), [1, 0])
        self.assertEqual(tables.select(info, tables.parse_selector("sample_type=pooled")), ["S2"])

    def test_process_and_qc_register_table_samples(self):
        for name in ("S1", "S2", "blank"):
            (self.root / f"{name}.mzML").touch()
        metadata = self.root / "metadata.csv"
        pd.DataFrame({"sample": ["S1", "S2", "blank"],
                      "sample_type": ["sample", "sample", "blank"]}).to_csv(metadata, index=False)
        out = self.root / "asari_output"
        out.mkdir()
        ft = pd.DataFrame({"id_number": [f"F{i}" for i in range(30)],
                           "mz": np.arange(30) + 100., "rtime": np.arange(30) + 10.,
                           "future_annotation": ["unknown"] * 30,
                           "S1": np.arange(30) + 100., "S2": np.arange(30) * 2 + 100.,
                           "blank": np.ones(30)})
        ft.to_csv(out / "preferred_Feature_table.tsv", sep="\t", index=False)
        args = cli.build_parser().parse_args(["process", "-p", str(self.root), "-i", str(self.root),
                                             "-m", "pos", "--metadata", str(metadata), "--single-run"])
        with patch("pcpfm.process.asari_runner.check_asari", return_value="1.18.0"), \
             patch("pcpfm.process._run_or_reuse", return_value={"output_dir": str(out)}):
            process._process(args, self.record)
        self.assertEqual(self.record.data["artifacts"]["feature_table"]["sample_columns"], ["S1", "S2"])
        self.assertEqual(self.record.data["artifacts"]["blank_feature_table"]["sample_columns"], ["blank"])
        args = cli.build_parser().parse_args(["qc", "-p", str(self.root)])
        with patch("pcpfm.report.write_report", return_value=str(self.root / "report.html")):
            qc._qc(args, self.record)
        for name in ("feature_table_filtered", "feature_table_analysis"):
            loaded, cols = tables.load_table(self.record, name)
            self.assertEqual(cols, ["S1", "S2"])
            self.assertEqual(loaded["future_annotation"].tolist(), ["unknown"] * 30)
        metrics = qc.sample_metrics

        def flag_second_sample(*args, **kwargs):
            smet, corr = metrics(*args, **kwargs)
            smet.at["S2", "flags"] = "low_feature_count"
            return smet, corr

        args.drop_flagged_samples = True
        with patch("pcpfm.qc.sample_metrics", side_effect=flag_second_sample), \
             patch("pcpfm.report.write_report", return_value=str(self.root / "report.html")):
            qc._qc(args, self.record)
        self.record.save()
        reloaded = ProjectRecord.load(str(self.root))
        for name in ("feature_table_filtered", "feature_table_analysis"):
            loaded, cols = tables.load_table(reloaded, name)
            self.assertEqual(cols, ["S1"])
            self.assertNotIn("S2", loaded.columns)
