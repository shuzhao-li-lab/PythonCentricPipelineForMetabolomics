"""Tests for ``pcpfm.annotate.gc_evidence``. Run with `python -m unittest tests.test_annotate_gc_evidence -v`."""

import json
import os
import tempfile
import unittest

import pandas as pd

from pcpfm.annotate.gc_evidence import (
    GC_EVIDENCE_COLUMNS,
    file_sha256,
    staged_dir,
    write_evidence,
    write_manifest,
)


class TestGcOutputs(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.out = os.path.join(tmp.name, "gc")
        os.mkdir(self.out)
        with open(os.path.join(self.out, "old.txt"), "w") as f:
            f.write("previous run")

    def test_success_replaces_outputs(self):
        with staged_dir(self.out) as staging:
            write_manifest({"workflow": "gc"}, os.path.join(staging, "manifest.json"))
        self.assertEqual(os.listdir(self.out), ["manifest.json"])
        self.assertEqual(os.listdir(os.path.dirname(self.out)), ["gc"])
        with open(os.path.join(self.out, "manifest.json")) as f:
            self.assertEqual(json.load(f), {"workflow": "gc"})

    def test_failure_keeps_previous_outputs(self):
        with self.assertRaises(RuntimeError), staged_dir(self.out) as staging:
            write_manifest({}, os.path.join(staging, "manifest.json"))
            raise RuntimeError("Asari failed")
        self.assertEqual(os.listdir(self.out), ["old.txt"])
        self.assertEqual(os.listdir(os.path.dirname(self.out)), ["gc"])

    def test_sha256_follows_content(self):
        path = os.path.join(self.out, "old.txt")
        before = file_sha256(path)
        with open(path, "w") as f:
            f.write("changed")
        self.assertNotEqual(file_sha256(path), before)


class TestGcEvidence(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = os.path.join(tmp.name, "annotations.tsv")

    def test_empty_is_header_only(self):
        self.assertEqual(write_evidence(pd.DataFrame(), self.path), 0)
        with open(self.path) as f:
            self.assertEqual(f.read().strip().split("\t"), GC_EVIDENCE_COLUMNS)

    def test_missing_columns_written_empty(self):
        rows = pd.DataFrame({"level": ["GC"], "feature_id": ["F1"]})
        self.assertEqual(write_evidence(rows, self.path), 1)
        out = pd.read_csv(self.path, sep="\t")
        self.assertEqual(list(out.columns), GC_EVIDENCE_COLUMNS)
        self.assertTrue(out["ri_error"].isna().all())

    def test_unknown_column_rejected(self):
        with self.assertRaises(ValueError):
            write_evidence(pd.DataFrame({"rt_error_s": [1.0]}), self.path)
