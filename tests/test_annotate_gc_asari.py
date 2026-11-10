"""Tests for ``pcpfm.annotate.gc_asari``. Run with `python -m unittest tests.test_annotate_gc_asari -v`."""

import json
import os
import tempfile
import unittest
from unittest import mock

import pandas as pd

from pcpfm.annotate.gc_asari import (
    denovo_outputs,
    gc_evidence_rows,
    run_asari_gc,
    write_asari_input,
)
from pcpfm.annotate.gc_evidence import GC_EVIDENCE_COLUMNS
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


# two rows of Asari's Features table: one empCpd matched through two features
FEATURES_TSV = (
    "feature\tmz\trtime\tRI\tempCpd\tquant_ion\tscore_cosine\tscore_entropy\tname\t"
    "inchikey\tlib_RI\tdelta_RI\tpeak_area\tcSelectivity\tpeak_shape\tsnr\tdetection_counts\n"
    "F75\t42.04\t414.29\t1188.0\tempCpd_0002\tF107\t0.397\t0.587\t3'-Hydroxy\"butanone\t"
    "NBEGXSQMVJTIAR-UHFFFAOYSA-N\t1178.0\t10.046\t1e8\t0.65\t0.91\t70.0\t4\n"
    "F1470\t71.05\t401.63\t1155.1\tempCpd_0002\tF107\t0.397\t0.587\tHydroxybutanone\t"
    "NBEGXSQMVJTIAR-UHFFFAOYSA-N\t1178.0\t-22.886\t2e8\t1.0\t0.94\t675.0\t4\n"
)


class TestEvidenceRows(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = os.path.join(tmp.name, "Features_gc_annotation.tsv")
        with open(self.path, "w") as f:
            f.write(FEATURES_TSV)

    def test_one_row_per_matched_feature(self):
        rows, dropped = gc_evidence_rows(self.path, "lib.msp", (300.0, 600.0))
        self.assertEqual(dropped, 0)
        self.assertEqual(list(rows.columns), GC_EVIDENCE_COLUMNS)
        self.assertEqual(list(rows["feature_id"]), ["F75", "F1470"])
        first = rows.iloc[0]
        self.assertEqual(first["empcpd_id"], "empCpd_0002")
        self.assertEqual(first["reference_id"], "NBEGXSQMVJTIAR-UHFFFAOYSA-N")
        self.assertEqual(first["reference_source"], "lib.msp")
        self.assertEqual(first["reference_name"], "3'-Hydroxy\"butanone")
        self.assertEqual(
            (first["level"], first["evidence_type"]), ("GC", "ei_ri_library")
        )
        self.assertEqual(list(rows["ri_error"]), [10.046, -22.886])
        self.assertEqual(first["score_entropy"], 0.587)

    def test_features_outside_kovats_range_dropped(self):
        # F75 at 414.29 s, F1470 at 401.63 s; the range bounds are inclusive
        rows, dropped = gc_evidence_rows(self.path, "lib.msp", (401.63, 410.0))
        self.assertEqual((list(rows["feature_id"]), dropped), (["F1470"], 1))
        self.assertEqual(list(rows.index), [0])

    def test_zero_hits_give_header_only_evidence(self):
        with open(self.path, "w") as f:
            f.write(FEATURES_TSV.splitlines()[0] + "\n")
        rows, dropped = gc_evidence_rows(self.path, "lib.msp", (300.0, 600.0))
        self.assertEqual((len(rows), dropped), (0, 0))
        self.assertEqual(list(rows.columns), GC_EVIDENCE_COLUMNS)


class TestDenovoOutputs(unittest.TestCase):
    def test_counts_denovo_pseudospectra(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(denovo_outputs(d), {})
            with open(
                os.path.join(d, "gc_annotation_denovo_pseudospectra.json"), "w"
            ) as f:
                json.dump([{"id": "ps1"}, {"id": "ps2"}], f)
            self.assertEqual(
                denovo_outputs(d),
                {
                    "json": "gc_annotation_denovo_pseudospectra.json",
                    "msp": "gc_annotation_denovo_pseudospectra.msp",
                    "n_pseudospectra": 2,
                },
            )


class TestRealAsari(unittest.TestCase):
    """Unmocked: real Asari matches three co-eluting, correlated fragments to a library entry."""

    PEAKS = ((73.0469, 999.0), (147.0657, 400.0), (191.0890, 250.0))

    def test_synthetic_compound_is_annotated(self):
        table = gc_table(n=3, samples=("S1", "S2", "S3"))
        table["mz"] = [mz for mz, _ in self.PEAKS]
        table["rtime"] = [450.0, 450.1, 450.2]  # RI 1250 with the Kovats file below
        for k, s in enumerate(("S1", "S2", "S3")):
            table[s] = [i * 1000 * 2**k for _, i in self.PEAKS]
        with tempfile.TemporaryDirectory() as d:
            kovats = os.path.join(d, "k.tsv")
            with open(kovats, "w") as f:
                f.write(
                    "RI\tRT\n"
                    + "".join(f"{1000 + 100 * j}\t{5 + j}\n" for j in range(6))
                )
            db = os.path.join(d, "lib.msp")
            with open(db, "w") as f:
                f.write("Name: TestCpd\nInChIKey: TEST-KEY\nRETENTIONTIME: 1250\n")
                f.write(f"ExactMass: 191.089\nNum Peaks: {len(self.PEAKS)}\n")
                f.write("".join(f"{mz} {i}\n" for mz, i in self.PEAKS) + "\n")
            result = run_asari_gc(table, ["S1", "S2", "S3"], kovats, db, PARAMS, d)
            rows, dropped = gc_evidence_rows(
                result["features_tsv"], "lib.msp", (300, 600)
            )
        self.assertEqual(sorted(rows["feature_id"]), ["F0", "F1", "F2"])
        self.assertEqual(set(rows["reference_id"]), {"TEST-KEY"})
        self.assertEqual(dropped, 0)
