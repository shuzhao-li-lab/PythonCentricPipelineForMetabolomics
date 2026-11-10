"""Tests for ``pcpfm.annotate.gc_evidence``. Run with `python -m unittest tests.test_annotate_gc_evidence -v`."""

import json
import os
import tempfile
import unittest

from pcpfm.annotate.gc_evidence import file_sha256, staged_dir, write_manifest


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
