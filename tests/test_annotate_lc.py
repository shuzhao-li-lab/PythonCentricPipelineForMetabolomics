"""Tests for ``pcpfm annotate lc``. Run with `python -m unittest tests.test_annotate_lc -v`."""

import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

from pcpfm import cli
from pcpfm.annotate.lc_evidence import (
    EVIDENCE_COLUMNS,
    empty_evidence,
    file_sha256,
    write_evidence,
)
from pcpfm.project import FAILED, SKIPPED, ProjectRecord


def make_project(d, workflow="LC", mode="pos"):
    """A project with a 5-feature, 3-sample feature_table and a recorded sample run."""
    record = ProjectRecord.load_or_create(d)
    path = os.path.join(record.subdir("tables"), "feature_table.tsv")
    pd.DataFrame(
        {
            "id_number": [f"F{i}" for i in range(5)],
            "mz": np.linspace(100, 500, 5),
            "rtime": np.linspace(30, 300, 5),
            "S1": 1.0,
            "S2": 2.0,
            "S3": 3.0,
        }
    ).to_csv(path, sep="\t", index=False)
    record.add_artifact(
        "feature_table", path, "test input", "process", sample_columns=["S1", "S2", "S3"]
    )
    record.data["runs"]["samples"] = {"workflow": workflow, "mode": mode}
    record.save()
    return record


def annotate(argv):
    """Run `pcpfm annotate ...` quietly; return what argparse wrote to stderr."""
    err = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
        try:
            cli.main(["annotate"] + argv)
        except SystemExit as e:
            raise SystemExit(err.getvalue()) from e
    return err.getvalue()


class TestEvidence(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = os.path.join(tmp.name, "annotations.tsv")

    def test_empty_is_header_only(self):
        self.assertEqual(write_evidence(empty_evidence(), self.path), 0)
        with open(self.path) as f:
            self.assertEqual(f.read().strip().split("\t"), EVIDENCE_COLUMNS)

    def test_missing_columns_written_empty_in_fixed_order(self):
        rows = pd.DataFrame({"level": ["4", "4"], "feature_id": ["F1", "F1"]})
        self.assertEqual(write_evidence(rows, self.path), 2)
        out = pd.read_csv(self.path, sep="\t", dtype=str)
        self.assertEqual(list(out.columns), EVIDENCE_COLUMNS)
        self.assertEqual(list(out["feature_id"]), ["F1", "F1"])
        self.assertTrue(out["score"].isna().all())

    def test_unknown_column_rejected(self):
        with self.assertRaises(ValueError):
            write_evidence(pd.DataFrame({"typo_score": [1.0]}), self.path)


class TestAnnotateLcCommand(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name

    def test_workflow_subcommand_required(self):
        make_project(self.d)
        with self.assertRaises(SystemExit) as cm:
            annotate([])
        self.assertIn("required: workflow", str(cm.exception))

    def test_foreign_workflow_flag_rejected(self):
        make_project(self.d)
        with self.assertRaises(SystemExit) as cm:
            annotate(["lc", "-p", self.d, "--ri-calibrants", "alkanes.tsv"])
        self.assertIn("unrecognized arguments: --ri-calibrants", str(cm.exception))

    def test_help_does_not_import_annotation(self):
        code = "import sys, pcpfm.cli; pcpfm.cli.build_parser(); print('pcpfm.annotate' in sys.modules)"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
        self.assertEqual(out.stdout.strip(), "False")

    def test_gc_source_fails_before_outputs(self):
        make_project(self.d, workflow="GC")
        with self.assertRaises(ValueError):
            annotate(["lc", "-p", self.d])
        record = ProjectRecord.load(self.d)
        self.assertEqual(record.step("annotate")["status"], FAILED)
        self.assertIn("GC", record.step("annotate")["message"])
        self.assertFalse(record.has_artifact("annotations_lc"))
        self.assertFalse(os.path.exists(os.path.join(self.d, "annotations")))
        self.assertEqual(record.data["history"][-1]["status"], "failed")

    def test_unknown_source_workflow_warns(self):
        make_project(self.d, workflow=None)
        annotate(["lc", "-p", self.d])
        record = ProjectRecord.load(self.d)
        self.assertIn("not recorded", record.step("annotate")["message"])
        with open(record.abspath(record.data["artifacts"]["annotations_lc"]["manifest"])) as f:
            manifest = json.load(f)
        self.assertEqual(len(manifest["warnings"]), 1)

    def test_outputs_and_record(self):
        make_project(self.d)
        annotate(["lc", "-p", self.d])
        record = ProjectRecord.load(self.d)
        step = record.step("annotate")
        # no level is implemented yet, so nothing ran
        self.assertEqual(step["status"], SKIPPED)
        self.assertEqual(step["params"]["workflow"], "lc")
        art = record.data["artifacts"]["annotations_lc"]
        evidence = pd.read_csv(record.abspath(art["path"]), sep="\t")
        self.assertEqual((len(evidence), list(evidence.columns)), (0, EVIDENCE_COLUMNS))
        with open(record.abspath(art["manifest"])) as f:
            manifest = json.load(f)
        source_sha = file_sha256(record.artifact_path("feature_table"))
        self.assertEqual(manifest["source"]["sha256"], source_sha)
        self.assertEqual(art["source_sha256"], source_sha)
        self.assertEqual(manifest["source"]["sample_columns"], ["S1", "S2", "S3"])
        self.assertEqual(manifest["source_context"], {"workflow": "LC", "mode": "pos"})
        self.assertEqual(step["summary"]["workflows"]["lc"]["run_id"], art["run_id"])

    def test_other_workflows_kept_in_summary(self):
        record = make_project(self.d)
        gc = {"run_id": "gc1", "levels": {"EI": {"status": "completed", "n_hits": 3}}}
        record.start_step("annotate")
        record.finish_step("annotate", summary={"workflows": {"gc": gc}})
        annotate(["lc", "-p", self.d])
        workflows = ProjectRecord.load(self.d).step("annotate")["summary"]["workflows"]
        self.assertEqual(workflows["gc"], gc)
        self.assertIn("lc", workflows)

    def test_failed_rerun_keeps_previous_success(self):
        make_project(self.d)
        annotate(["lc", "-p", self.d])
        first = ProjectRecord.load(self.d)
        gc = {"run_id": "gc1", "levels": {}}
        first.data["steps"]["annotate"]["summary"]["workflows"]["gc"] = gc
        art = first.data["artifacts"]["annotations_lc"]
        first.data["runs"]["samples"]["workflow"] = "GC"
        first.save()
        with self.assertRaises(ValueError):
            annotate(["lc", "-p", self.d])
        record = ProjectRecord.load(self.d)
        self.assertEqual(record.data["artifacts"]["annotations_lc"], art)
        workflows = record.step("annotate")["summary"]["workflows"]
        self.assertEqual(workflows["lc"]["run_id"], art["run_id"])
        self.assertEqual(workflows["gc"], gc)

    def test_status_shows_workflow_levels(self):
        make_project(self.d)
        annotate(["lc", "-p", self.d])
        record = ProjectRecord.load(self.d)
        record.data["steps"]["annotate"]["summary"]["workflows"]["lc"]["levels"]["L4"] = {
            "status": "completed",
            "n_features": 4,
            "n_hits": 9,
        }
        text = "\n".join(record.status_lines())
        self.assertIn("    lc (run ", text)
        self.assertIn("L4 4 features, 9 hits", text)
        self.assertIn("L1b skipped", text)


if __name__ == "__main__":
    unittest.main()
