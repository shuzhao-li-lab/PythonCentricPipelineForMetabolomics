"""Tests for ``pcpfm annotate lc``. Run with `python -m unittest tests.test_annotate_lc -v`."""

import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import types
import unittest
import warnings
from unittest import mock

import numpy as np
import pandas as pd

from pcpfm import cli
from pcpfm.annotate import lc
from pcpfm.annotate.lc_evidence import (
    EVIDENCE_COLUMNS,
    empty_evidence,
    file_sha256,
    write_evidence,
)
from pcpfm.annotate.lc_preannotation import (
    ASARI_OUTPUT,
    L4_SOURCE,
    empcpd_membership,
    l4_evidence,
    mz_correction,
    run_asari_preannotation,
    write_asari_input,
)
from pcpfm.project import COMPLETED, FAILED, ProjectRecord
from pcpfm.tables import FEATURE_COLUMNS


# Asari's three empCpd kinds; F1 belongs to the khipu empCpd and to the singleton
EMPCPDS = {
    "kp1_145.1577": {
        "interim_id": "kp1_145.1577",
        "neutral_formula_mass": 145.157897623 * (1 + 2e-6),
        "neutral_formula": None,
        "list_matches": [["C7H19N3_145.157898", "", 1], ["ISO_145.158", "", 1]],
        "MS1_pseudo_Spectra": [
            {"id_number": "F0", "ion_relation": "M0,M+H+"},
            {"id_number": "F1", "ion_relation": "13C/12C,M+H+"},
        ],
    },
    "kp2_200.0": {
        "interim_id": "kp2_200.0",
        "neutral_formula_mass": 200.0,
        "neutral_formula": None,
        "list_matches": [],
        "MS1_pseudo_Spectra": [{"id_number": "F2", "ion_relation": "M0,M+H+"}],
    },
    "7": {
        "interim_id": "_singleton_7",
        "neutral_formula_mass": 126.03,
        "neutral_formula": "C6H6O3",
        "MS1_pseudo_Spectra": [{"id_number": "F3"}, {"id_number": "F1", "ion_relation": "M+Na+"}],
    },
    "100003": {
        "interim_id": 100003,
        "neutral_formula_mass": "",
        "neutral_formula": "",
        "MS1_pseudo_Spectra": [{"id_number": "F4"}],
    },
}
COMPOUNDS = {
    "C7H19N3_145.157898": {
        "neutral_formula": "C7H19N3",
        "neutral_formula_mass": 145.157897623,
        "compounds": [{"primary_id": "HMDB0001257", "name": "Spermidine"}],
    },
    "ISO_145.158": {
        "neutral_formula": "ISO",
        "neutral_formula_mass": 145.158,
        "compounds": [{"primary_id": "H1", "name": "A"}, {"primary_id": "H2", "name": "B"}],
    },
}


def make_project(d, workflow="LC", mode="pos"):
    """A project with a 5-feature, 3-sample feature_table and a recorded sample run."""
    record = ProjectRecord.load_or_create(d)
    path = os.path.join(record.subdir("tables"), "feature_table.tsv")
    asari_table(5, ["S1", "S2", "S3"]).to_csv(path, sep="\t", index=False)
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
        patch = fake_asari()
        patch.start()
        self.addCleanup(patch.stop)

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
        self.assertEqual(FakeAnnotation.calls, [])

    def test_missing_mode_fails(self):
        make_project(self.d, mode=None)
        with self.assertRaisesRegex(ValueError, "pos' or 'neg"):
            annotate(["lc", "-p", self.d])
        record = ProjectRecord.load(self.d)
        self.assertEqual(record.step("annotate")["status"], FAILED)
        self.assertFalse(record.has_artifact("annotations_lc"))
        self.assertEqual(FakeAnnotation.calls, [])

    def test_preannotation_flags_and_mode_reach_asari(self):
        make_project(self.d, mode="neg")
        annotate(["lc", "-p", self.d, "--preannot-ppm", "7", "--preannot-rt", "3.5"])
        params, table = FakeAnnotation.calls[0]
        self.assertEqual(
            {k: params[k] for k in ("mode", "mz_tolerance_ppm", "khipu_rtime_tolerance")},
            {"mode": "neg", "mz_tolerance_ppm": 7.0, "khipu_rtime_tolerance": 3.5},
        )
        self.assertEqual(list(table.columns), FEATURE_COLUMNS + ["S1", "S2", "S3"])

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
        self.assertEqual(step["status"], COMPLETED)
        self.assertEqual(step["params"]["workflow"], "lc")
        self.assertEqual((step["params"]["preannot_ppm"], step["params"]["preannot_rt"]), (5, 2))
        levels = step["summary"]["workflows"]["lc"]["levels"]
        self.assertEqual(
            levels["grouping"],
            {"status": COMPLETED, "n_empcpds": 2, "n_singletons": 1, "n_orphans": 1},
        )
        self.assertEqual(
            levels["L4"], {"status": COMPLETED, "n_empcpds": 1, "n_features": 2, "n_rows": 4}
        )
        art = record.data["artifacts"]["annotations_lc"]
        out_dir = os.path.join(self.d, "annotations", "lc")
        self.assertEqual(record.abspath(art["path"]), os.path.join(out_dir, "annotations.tsv"))
        self.assertTrue(os.path.exists(os.path.join(out_dir, ASARI_OUTPUT)))
        evidence = pd.read_csv(record.abspath(art["path"]), sep="\t")
        self.assertEqual((len(evidence), list(evidence.columns)), (4, EVIDENCE_COLUMNS))
        self.assertEqual(art["n_rows"], 4)
        with open(record.abspath(art["manifest"])) as f:
            manifest = json.load(f)
        source_sha = file_sha256(record.artifact_path("feature_table"))
        self.assertEqual(manifest["source"]["sha256"], source_sha)
        self.assertEqual(art["source_sha256"], source_sha)
        self.assertEqual(manifest["source"]["sample_columns"], ["S1", "S2", "S3"])
        self.assertEqual(manifest["source_context"], {"workflow": "LC", "mode": "pos"})
        self.assertEqual(
            manifest["preannotation"],
            {
                "asari_version": "9.9.9",
                "asari_params": {"mode": "pos", "mz_tolerance_ppm": 5, "khipu_rtime_tolerance": 2},
                "mz_corrected_by_division": None,
            },
        )
        self.assertEqual(manifest["outputs"]["empirical_compounds"], ASARI_OUTPUT)
        self.assertEqual(manifest["levels"], levels)

    def test_rerun_replaces_outputs(self):
        make_project(self.d)
        annotate(["lc", "-p", self.d])
        out_dir = os.path.join(self.d, "annotations", "lc")
        with open(os.path.join(out_dir, "leftover.txt"), "w") as f:
            f.write("x")
        annotate(["lc", "-p", self.d, "--preannot-ppm", "7"])
        self.assertEqual(os.listdir(os.path.join(self.d, "annotations")), ["lc"])
        self.assertNotIn("leftover.txt", os.listdir(out_dir))
        with open(os.path.join(out_dir, "manifest.json")) as f:
            self.assertEqual(json.load(f)["params"]["preannot_ppm"], 7)

    def test_asari_error_keeps_previous_run_and_names_log(self):
        make_project(self.d)
        annotate(["lc", "-p", self.d])
        art = ProjectRecord.load(self.d).data["artifacts"]["annotations_lc"]
        with mock.patch.object(FakeAnnotation, "error", KeyError("mz")):
            with self.assertRaisesRegex(RuntimeError, "KeyError.*fake asari console output"):
                annotate(["lc", "-p", self.d])
        record = ProjectRecord.load(self.d)
        self.assertEqual(record.step("annotate")["status"], FAILED)
        self.assertIn("Asari pre-annotation failed", record.step("annotate")["message"])
        self.assertEqual(record.data["artifacts"]["annotations_lc"], art)
        # the previous outputs survive and no staging directory is left behind
        self.assertEqual(os.listdir(os.path.join(self.d, "annotations")), ["lc"])
        evidence = pd.read_csv(record.abspath(art["path"]), sep="\t")
        self.assertEqual(len(evidence), 4)

    def test_failed_optional_step_does_not_stop_the_run(self):
        make_project(self.d)
        later = mock.Mock(return_value=(empty_evidence(), {"L2": {"status": COMPLETED, "n_rows": 0}}))
        seen = []

        def failing(run):
            seen.append(len(run.membership))
            raise ValueError("bad standards file")

        plan = lc._plan
        steps = lambda args: plan(args) + [
            lc.Step("L1b", failing, required=False),
            lc.Step("L2", later, required=False),
        ]
        with mock.patch.object(lc, "_plan", steps):
            annotate(["lc", "-p", self.d])
        # later steps see what pre-annotation produced
        self.assertEqual(seen, [6])
        record = ProjectRecord.load(self.d)
        step = record.step("annotate")
        self.assertEqual(step["status"], COMPLETED)
        self.assertIn("L1b failed: bad standards file", step["message"])
        levels = step["summary"]["workflows"]["lc"]["levels"]
        self.assertEqual(levels["L1b"], {"status": FAILED, "reason": "bad standards file"})
        self.assertEqual(levels["L2"]["status"], COMPLETED)
        self.assertEqual(levels["L4"]["n_rows"], 4)
        self.assertEqual(record.data["artifacts"]["annotations_lc"]["n_rows"], 4)

    def test_optional_step_with_bad_columns_fails_alone(self):
        make_project(self.d)
        bad = mock.Mock(return_value=(pd.DataFrame({"typo": [1]}), {}))
        plan = lc._plan
        with mock.patch.object(lc, "_plan", lambda args: plan(args) + [lc.Step("L2", bad, False)]):
            annotate(["lc", "-p", self.d])
        step = ProjectRecord.load(self.d).step("annotate")
        self.assertEqual(step["status"], COMPLETED)
        self.assertEqual(step["summary"]["workflows"]["lc"]["levels"]["L2"]["status"], FAILED)

    def test_failed_preannotation_stops_later_steps(self):
        make_project(self.d)
        later = mock.Mock()
        plan = lc._plan
        steps = lambda args: plan(args) + [lc.Step("L1b", later, required=False)]
        with mock.patch.object(lc, "_plan", steps), mock.patch.object(
            FakeAnnotation, "error", KeyError("mz")
        ):
            with self.assertRaisesRegex(RuntimeError, "Asari pre-annotation failed"):
                annotate(["lc", "-p", self.d])
        later.assert_not_called()
        record = ProjectRecord.load(self.d)
        self.assertEqual(record.step("annotate")["status"], FAILED)
        self.assertIn("Asari pre-annotation failed", record.step("annotate")["message"])
        self.assertFalse(record.has_artifact("annotations_lc"))

    def test_output_dir_not_private(self):
        make_project(self.d)
        annotate(["lc", "-p", self.d])
        umask = os.umask(0)
        os.umask(umask)
        mode = os.stat(os.path.join(self.d, "annotations", "lc")).st_mode & 0o777
        self.assertEqual(mode, 0o777 & ~umask)

    def test_failed_swap_restores_previous_outputs(self):
        make_project(self.d)
        annotate(["lc", "-p", self.d])
        out_dir = os.path.join(self.d, "annotations", "lc")
        before = sorted(os.listdir(out_dir))
        with mock.patch("pcpfm.annotate.lc_evidence.os.replace", side_effect=OSError("disk full")):
            with self.assertRaises(OSError):
                annotate(["lc", "-p", self.d])
        self.assertEqual(os.listdir(os.path.join(self.d, "annotations")), ["lc"])
        self.assertEqual(sorted(os.listdir(out_dir)), before)

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
        self.assertEqual(workflows["lc"]["levels"]["L4"]["n_rows"], 4)
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
        self.assertIn("    lc: grouping 2 empcpds, 1 singletons, 1 orphans", text)
        self.assertIn("L4 4 features, 9 hits", text)
        self.assertIn("L1b skipped", text)


def asari_table(n=4, samples=("S1", "S2")):
    """A pcpfm feature table with all Asari feature columns, a blank-filter column and samples."""
    t = pd.DataFrame({"id_number": [f"F{i}" for i in range(n)]})
    for c in FEATURE_COLUMNS[1:]:
        t[c] = np.arange(n, dtype=float)
    t["blank_flag"] = False
    for i, c in enumerate(samples, start=1):
        t[c] = float(i)
    return t


class FakeAnnotation:
    """Stands in for asari.annotate.LCMS_Annotation; records what it was given."""

    calls = []
    writes_output = True
    error = None

    def __init__(self, parameters):
        self.parameters = parameters

    def annotate_user_featuretable(self, infile):
        table = pd.read_csv(infile, sep="\t", dtype={"id_number": str})
        FakeAnnotation.calls.append((dict(self.parameters), table))
        FakeAnnotation.infile = infile
        print("fake asari console output")
        if self.error:
            raise self.error
        self.KCD = types.SimpleNamespace(mass_indexed_compounds=COMPOUNDS)
        if self.writes_output:
            with open(os.path.join(self.parameters["outdir"], ASARI_OUTPUT), "w") as f:
                json.dump(EMPCPDS, f)


def fake_asari():
    """Patch sys.modules so the adapter imports FakeAnnotation instead of Asari."""
    FakeAnnotation.calls = []
    pkg = types.ModuleType("asari")
    pkg.__version__ = "9.9.9"
    mod = types.ModuleType("asari.annotate")
    mod.LCMS_Annotation = FakeAnnotation
    pkg.annotate = mod
    return mock.patch.dict(sys.modules, {"asari": pkg, "asari.annotate": mod})


class TestAsariInput(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.path = os.path.join(tmp.name, "in.tsv")

    def test_asari_column_order_without_other_columns(self):
        t = asari_table()
        write_asari_input(t[["S2", "blank_flag"] + FEATURE_COLUMNS[::-1] + ["S1"]], ["S1", "S2"], self.path)
        with open(self.path) as f:
            header, first = f.readline().rstrip("\n").split("\t"), f.readline().split("\t")
        self.assertEqual(header, FEATURE_COLUMNS + ["S1", "S2"])
        # Asari parses these two with int()
        self.assertEqual(first[FEATURE_COLUMNS.index("parent_masstrack_id")], "0")
        self.assertEqual(first[FEATURE_COLUMNS.index("detection_counts")], "0")

    def test_real_asari_reader_parses_output(self):
        try:
            from asari.tools.file_io import read_features_from_asari_table
        except ImportError:
            self.skipTest("asari is not installed")
        t = asari_table()
        t["snr"] = [3.5, 4.0, 5.0, 6.0]
        write_asari_input(t, ["S1", "S2"], self.path)
        with open(self.path) as f, contextlib.redirect_stdout(io.StringIO()):
            n_samples, peaks = read_features_from_asari_table(f.read())
        self.assertEqual(n_samples, 2)
        self.assertEqual([p["id_number"] for p in peaks], ["F0", "F1", "F2", "F3"])
        self.assertEqual(peaks[2]["mz"], 2.0)
        self.assertEqual(peaks[1]["snr"], 4.0)
        self.assertEqual(peaks[3]["parent_masstrack_id"], 3)

    def test_infinite_value(self):
        t = asari_table()
        t.loc[1, "mz"] = np.inf
        with self.assertRaisesRegex(ValueError, "mz"):
            write_asari_input(t, ["S1", "S2"], self.path)

    def test_sample_named_like_feature_column(self):
        t = asari_table().rename(columns={"S1": "snr_dup"})
        t["snr_dup"] = 1.0
        with self.assertRaisesRegex(ValueError, "snr"):
            write_asari_input(t, ["snr", "S2"], self.path)

    def test_missing_feature_column(self):
        with self.assertRaisesRegex(ValueError, "cSelectivity"):
            write_asari_input(asari_table().drop(columns="cSelectivity"), ["S1", "S2"], self.path)

    def test_missing_value(self):
        t = asari_table()
        t.loc[1, "snr"] = np.nan
        with self.assertRaisesRegex(ValueError, "snr"):
            write_asari_input(t, ["S1", "S2"], self.path)

    def test_non_integer_track_id(self):
        t = asari_table()
        t.loc[1, "parent_masstrack_id"] = 1.5
        with self.assertRaisesRegex(ValueError, "parent_masstrack_id"):
            write_asari_input(t, ["S1", "S2"], self.path)


class TestAsariPreannotation(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name

    def test_params_per_mode(self):
        for mode in ("pos", "neg"):
            with fake_asari():
                out = run_asari_preannotation(asari_table(), ["S1", "S2"], mode, 7, 3, self.d)
            params, table = FakeAnnotation.calls[0]
            self.assertEqual(
                params,
                {"mode": mode, "mz_tolerance_ppm": 7, "khipu_rtime_tolerance": 3, "outdir": self.d},
            )
            self.assertEqual(list(table["id_number"]), ["F0", "F1", "F2", "F3"])
            self.assertEqual(out["asari_version"], "9.9.9")
            self.assertEqual(out["empcpds_json"], os.path.join(self.d, ASARI_OUTPUT))
            with open(out["log"]) as f:
                self.assertIn("fake asari console output", f.read())

    def test_asari_input_copy_removed(self):
        with fake_asari():
            run_asari_preannotation(asari_table(), ["S1", "S2"], "pos", 5, 2, self.d)
        self.assertFalse(os.path.exists(FakeAnnotation.infile))
        self.assertNotIn("asari_input.tsv", os.listdir(self.d))

    def test_asari_input_copy_removed_on_failure(self):
        with fake_asari(), mock.patch.object(FakeAnnotation, "writes_output", False):
            with self.assertRaises(RuntimeError):
                run_asari_preannotation(asari_table(), ["S1", "S2"], "pos", 5, 2, self.d)
        self.assertFalse(os.path.exists(FakeAnnotation.infile))

    def test_unknown_mode_fails_before_asari(self):
        with fake_asari(), self.assertRaisesRegex(ValueError, "pos' or 'neg"):
            run_asari_preannotation(asari_table(), ["S1", "S2"], None, 5, 2, self.d)
        self.assertEqual(FakeAnnotation.calls, [])

    def test_missing_output(self):
        # a leftover JSON from an earlier call must not count as output
        with open(os.path.join(self.d, ASARI_OUTPUT), "w") as f:
            f.write("{}")
        with fake_asari(), mock.patch.object(FakeAnnotation, "writes_output", False):
            with self.assertRaises(RuntimeError):
                run_asari_preannotation(asari_table(), ["S1", "S2"], "pos", 5, 2, self.d)


class TestEmpcpdParsing(unittest.TestCase):
    def test_membership_kinds_and_shared_feature(self):
        m = empcpd_membership(EMPCPDS)
        self.assertEqual(len(m), 6)
        kinds = dict(zip(m["empcpd_id"], m["kind"]))
        self.assertEqual(
            kinds,
            {"kp1_145.1577": "khipu", "kp2_200.0": "khipu", "_singleton_7": "singleton", "100003": "orphan"},
        )
        self.assertEqual(sorted(m.loc[m["feature_id"] == "F1", "empcpd_id"]), ["_singleton_7", "kp1_145.1577"])

    def test_l4_one_row_per_feature_and_formula(self):
        rows = l4_evidence(EMPCPDS, empcpd_membership(EMPCPDS), COMPOUNDS)
        self.assertEqual(len(rows), 4)
        self.assertFalse(rows.duplicated(["feature_id", "empcpd_id", "reference_id"]).any())
        # the singleton membership of F1 adds no L4 row
        self.assertEqual(set(rows["empcpd_id"]), {"kp1_145.1577"})
        spd = rows[(rows["feature_id"] == "F1") & (rows["neutral_formula"] == "C7H19N3")].iloc[0]
        self.assertEqual((spd["reference_id"], spd["reference_name"]), ("HMDB0001257", "Spermidine"))
        self.assertEqual(spd["ion_relation"], "13C/12C,M+H+")
        self.assertAlmostEqual(spd["mz_error_ppm"], 2.0, places=6)
        iso = rows[rows["neutral_formula"] == "ISO"].iloc[0]
        self.assertEqual((iso["reference_id"], iso["reference_name"]), ("H1;H2", "A;B"))
        self.assertEqual(set(rows["level"]), {"4"})
        self.assertEqual(set(rows["reference_source"]), {L4_SOURCE})

    def test_l4_without_matches(self):
        empcpds = {"kp2_200.0": EMPCPDS["kp2_200.0"]}
        rows = l4_evidence(empcpds, empcpd_membership(empcpds), COMPOUNDS)
        self.assertEqual((len(rows), list(rows.columns)), (0, EVIDENCE_COLUMNS))

    def test_l4_match_missing_from_database(self):
        with self.assertRaisesRegex(ValueError, "ISO_145.158"):
            l4_evidence(EMPCPDS, empcpd_membership(EMPCPDS), {"C7H19N3_145.157898": COMPOUNDS["C7H19N3_145.157898"]})

    def test_mz_correction(self):
        self.assertIsNone(mz_correction(EMPCPDS))
        corrected = {"kp": {"MS1_pseudo_Spectra": [{"id_number": "F0", "mz_corrected_by_division": 1.000003}]}}
        self.assertEqual(mz_correction(corrected), 1.000003)


class TestRealAsari(unittest.TestCase):
    def test_spermidine_ions_grouped_and_matched(self):
        try:
            import asari.annotate  # noqa: F401
        except ImportError:
            self.skipTest("asari is not installed")
        # spermidine M+H+, its 13C isotope and M+Na+ at one RT, plus an unrelated feature
        t = asari_table()
        t["mz"] = [146.165174, 147.168529, 168.147116, 250.1234]
        t["rtime"] = [100.0, 100.5, 100.2, 300.0]
        t["peak_area"] = t["S1"] = t["S2"] = [1e7, 1e6, 3e6, 1e5]
        # khipu writes khipu.log to the working directory
        with tempfile.TemporaryDirectory() as d, contextlib.chdir(d), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = run_asari_preannotation(t, ["S1", "S2"], "pos", 5, 2, d)
            with open(out["empcpds_json"]) as f:
                empcpds = json.load(f)
        membership = empcpd_membership(empcpds)
        rows = l4_evidence(empcpds, membership, out["compounds"])
        spd = rows[rows["reference_name"] == "Spermidine"]
        self.assertEqual(sorted(spd["feature_id"]), ["F0", "F1", "F2"])
        self.assertEqual(spd["empcpd_id"].nunique(), 1)
        self.assertLess(abs(spd["mz_error_ppm"].iloc[0]), 1)
        self.assertNotIn("F3", set(rows["feature_id"]))
        self.assertIn("F3", set(membership["feature_id"]))


if __name__ == "__main__":
    unittest.main()
