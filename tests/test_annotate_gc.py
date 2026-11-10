"""Tests for ``pcpfm annotate gc``. Run with `python -m unittest tests.test_annotate_gc -v`."""

import contextlib
import io
import json
import os
import tempfile
import unittest
from unittest import mock

import pandas as pd

from pcpfm import cli
from pcpfm.annotate import gc
from pcpfm.project import ProjectRecord
from tests.test_annotate_gc_asari import gc_table

REQUIRED = ["annotate", "gc", "-p", "proj", "--kovats", "k.tsv", "--db", "lib.msp"]


def parse(argv):
    """Parse quietly; on error return argparse's message as the SystemExit text."""
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        try:
            return cli.build_parser().parse_args(argv)
        except SystemExit as e:
            raise SystemExit(err.getvalue()) from e


class TestAnnotateGcFlags(unittest.TestCase):
    def test_defaults_follow_asari(self):
        args = parse(REQUIRED)
        self.assertEqual(args.table, "feature_table")
        self.assertEqual(args.denovo, "F")
        self.assertEqual(args.denovo_limit, 20000)
        self.assertEqual(args.ri_tolerance, 30)
        self.assertEqual(args.score_cutoff_cosine, 0.4)
        self.assertEqual(args.score_cutoff_entropy, 0.4)
        self.assertEqual(args.corr_cutoff, 0.6)

    def test_kovats_and_db_required(self):
        for flag in ("--kovats", "--db"):
            argv = list(REQUIRED)
            i = argv.index(flag)
            del argv[i : i + 2]
            with self.assertRaises(SystemExit) as cm:
                parse(argv)
            self.assertIn(flag, str(cm.exception))

    def test_denovo_accepts_only_t_or_f(self):
        self.assertEqual(parse(REQUIRED + ["--denovo", "T"]).denovo, "T")
        with self.assertRaises(SystemExit):
            parse(REQUIRED + ["--denovo", "True"])

    def test_denovo_limit_must_be_positive_integer(self):
        self.assertEqual(parse(REQUIRED + ["--denovo-limit", "50"]).denovo_limit, 50)
        for bad in ("0", "-5", "2.5"):
            with self.assertRaises(SystemExit):
                parse(REQUIRED + ["--denovo-limit", bad])

    def test_lc_flags_rejected(self):
        with self.assertRaises(SystemExit) as cm:
            parse(REQUIRED + ["--ms1-standards", "std.csv"])
        self.assertIn("unrecognized arguments: --ms1-standards", str(cm.exception))

    def test_not_implemented_yet(self):
        with self.assertRaises(SystemExit) as cm:
            cli.main(REQUIRED)
        self.assertIn("not implemented yet", str(cm.exception))


class TestSourceContext(unittest.TestCase):
    def record(self, d, workflow):
        record = ProjectRecord.load_or_create(d)
        record.data["runs"]["samples"] = {"workflow": workflow, "mode": "pos"}
        return record

    def test_requires_gc_workflow(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(
                gc.source_context(self.record(d, "GC")), {"workflow": "GC"}
            )
        for workflow, found in (("LC", "Asari LC workflow"), (None, "unrecorded")):
            with (
                tempfile.TemporaryDirectory() as d,
                self.assertRaisesRegex(ValueError, found),
            ):
                gc.source_context(self.record(d, workflow))


# Asari output for the fake run: F0 matched inside the Kovats range, F2 outside it
FAKE_FEATURES = (
    "feature\tmz\trtime\tRI\tempCpd\tquant_ion\tscore_cosine\tscore_entropy\tname\t"
    "inchikey\tlib_RI\tdelta_RI\n"
    "F0\t73.0\t400.0\t1150.0\tempCpd_0001\tF0\t0.8\t0.7\tA\tKEY-A\t1140.0\t10.0\n"
    "F2\t75.0\t900.0\t1900.0\tempCpd_0002\tF2\t0.8\t0.7\tB\tKEY-B\t1900.0\t0.0\n"
)


def fake_asari(infile, outdir, KovatsIndex, database_file, **kwargs):
    """Stand-in for Asari's GC annotation: write the fake Features table."""
    name = kwargs["project_name_handle"]
    with open(os.path.join(outdir, f"Features_{name}.tsv"), "w") as f:
        f.write(FAKE_FEATURES)


def make_gc_project(d):
    """A GC project with a 3-feature table, a Kovats file and a 1-entry library."""
    record = ProjectRecord.load_or_create(d)
    path = os.path.join(record.subdir("tables"), "feature_table.tsv")
    gc_table().to_csv(path, sep="\t", index=False)
    record.add_artifact(
        "feature_table", path, "test input", "process", sample_columns=["S1", "S2"]
    )
    record.data["runs"]["samples"] = {"workflow": "GC", "mode": "pos"}
    record.save()
    with open(os.path.join(d, "kovats.tsv"), "w") as f:
        f.write("RI\tRT\n1000\t5.0\n1500\t10.0\n")  # 300-600 s
    with open(os.path.join(d, "lib.msp"), "w") as f:
        f.write(
            "Name: A\nRETENTIONTIME: 1140\nExactMass: 100.0\nNum Peaks: 1\n73.0 999\n\n"
        )
    return record


def gc_args(d):
    """Parsed `annotate gc` options for the test project."""
    k, db = os.path.join(d, "kovats.tsv"), os.path.join(d, "lib.msp")
    return cli.build_parser().parse_args(
        ["annotate", "gc", "-p", d, "--kovats", k, "--db", db]
    )


class TestAnnotateGc(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name
        self.record = make_gc_project(self.d)
        patch = mock.patch("asari.annotate.annotate_gcms_full", side_effect=fake_asari)
        patch.start()
        self.addCleanup(patch.stop)

    def out(self, name):
        return os.path.join(self.d, "annotations", "gc", name)

    def test_writes_evidence_and_manifest(self):
        with contextlib.redirect_stdout(io.StringIO()):
            summary, warnings = gc._annotate(gc_args(self.d), self.record, {})
        rows = pd.read_csv(self.out("annotations.tsv"), sep="\t")
        self.assertEqual(list(rows["feature_id"]), ["F0"])
        self.assertEqual(
            (summary["n_evidence_rows"], summary["n_outside_kovats"]), (1, 1)
        )
        self.assertIn("outside the Kovats range", warnings[0])
        with open(self.out("manifest.json")) as f:
            manifest = json.load(f)
        self.assertEqual(manifest["source_context"], {"workflow": "GC"})
        self.assertEqual(manifest["counts"]["n_library_entries"], 1)
        self.assertEqual(manifest["references"]["kovats"]["rt_range_s"], [300.0, 600.0])
        self.assertEqual(manifest["outputs"]["denovo"], {})
        artifact = ProjectRecord.load(self.d).data["artifacts"]["annotations_gc"]
        self.assertEqual(artifact["manifest"], "annotations/gc/manifest.json")
        self.assertEqual(
            (artifact["n_rows"], artifact["source_table"]), (1, "feature_table")
        )
        self.assertEqual(artifact["source_sha256"], manifest["source"]["sha256"])
