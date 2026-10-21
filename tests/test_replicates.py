"""Replicate correlation, flagging and collapse (`pcpfm replicates`)."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from pcpfm import cli, replicates
from pcpfm.project import COMPLETED, FAILED, ProjectRecord

N = 300


class TestAssessGroup(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        self.x, self.y = rng.normal(size=N), rng.normal(size=N)

    def test_keep_replicates_with_one_agreeing_pair(self):
        # r(1,2) = r(1,3) ~ 0.71, r(2,3) ~ 0: every replicate agrees with Rep1
        L = pd.DataFrame({"a": self.x + self.y, "b": self.x, "c": self.y})
        res = replicates.assess_group(L, ["a", "b", "c"], 0.6, drop_flagged=True)
        self.assertEqual(res["flag"], "min < 0.6")
        self.assertEqual(res["kept"], ["a", "b", "c"])
        self.assertEqual(res["representative"], "a")
        self.assertEqual(
            set(res["pairs"]), {"Rep1 vs Rep2", "Rep1 vs Rep3", "Rep2 vs Rep3"}
        )

    def test_failed_replicate_constant(self):
        L = pd.DataFrame({"a": np.zeros(N), "b": self.x, "c": self.x + 0.1 * self.y})
        res = replicates.assess_group(L, ["a", "b", "c"], 0.7, drop_flagged=True)
        self.assertEqual(res["flag"], "min < 0.7")
        self.assertEqual(res["dropped"], ["a"])
        self.assertEqual(res["representative"], "b")
        self.assertTrue(np.isfinite(res["min"]))

    def test_no_agreeing_pair(self):
        L = pd.DataFrame({"a": self.x, "b": self.y})
        res = replicates.assess_group(L, ["a", "b"], 0.7, drop_flagged=True)
        self.assertEqual(res["flag"], "max < 0.7")
        self.assertEqual(res["status"], replicates.STATUS_EXCLUDED)
        self.assertEqual(res["kept"], [])
        kept = replicates.assess_group(L, ["a", "b"], 0.7, drop_flagged=False)
        self.assertEqual(kept["flag"], "max < 0.7")
        self.assertEqual(kept["kept"], ["a", "b"])

    def test_single_replicate(self):
        res = replicates.assess_group(pd.DataFrame({"a": self.x}), ["a"], 0.7, True)
        self.assertEqual((res["flag"], res["kept"]), ("", ["a"]))
        self.assertTrue(np.isnan(res["min"]))


class TestReplicatesCommand(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.record = ProjectRecord.load_or_create(str(self.root))

        rng = np.random.default_rng(0)
        base = rng.lognormal(10, 2, N)
        # "rep": base intensities with 10% noise; "noise": unrelated intensities
        kinds = {
            "A1": "rep",
            "A2": "rep",
            "A3": "rep",  # good triplicate
            "B1": "rep",
            "B2": "rep",
            "B3": "noise",  # B3 failed
            "C1": "noise",
            "C2": "noise",  # no agreement: excluded
            "D1": "noise",
            "D2": "rep",
            "D3": "rep",  # first injection failed
            "E1": "rep",  # single replicate
            "QC1": "rep",  # pooled QC, no replicate group
        }
        cols = {
            s: base * rng.lognormal(0, 0.1, N)
            if k == "rep"
            else rng.lognormal(10, 2, N)
            for s, k in kinds.items()
        }
        self.samples = list(cols)
        ft = pd.DataFrame(
            {"id_number": [f"F{i}" for i in range(N)], "mz": np.linspace(100, 900, N)}
        )
        ft["blank_flag"] = False
        ft = pd.concat([ft, pd.DataFrame(cols)], axis=1)
        self.ft = ft
        tables = self.root / "tables"
        tables.mkdir()
        ft.to_csv(tables / "feature_table.tsv", sep="\t", index=False)
        self.record.add_artifact(
            "feature_table",
            tables / "feature_table.tsv",
            "test",
            "process",
            sample_columns=self.samples,
        )
        groups = {"QC1": ""}
        groups.update({s: s[0] for s in self.samples if s != "QC1"})
        # B is injected in reverse table order
        order = {"B1": 12, "B2": 11, "B3": 10}
        self.meta = pd.DataFrame(
            {
                "sample": self.samples,
                "sample_type": ["sample"] * 12 + ["pooled"],
                "biosample": [groups[s] for s in self.samples],
                "injection_order": [
                    order.get(s, i) for i, s in enumerate(self.samples)
                ],
            }
        )
        self.meta.to_csv(self.root / "metadata.tsv", sep="\t", index=False)
        self.record.add_artifact(
            "metadata", self.root / "metadata.tsv", "test", "process"
        )
        self.record.save()

    def run_cli(self, *extra):
        cli.main(["replicates", "-p", str(self.root), "--rep-col", "biosample", *extra])
        return ProjectRecord.load(str(self.root))

    def test_collapse_and_report(self):
        record = self.run_cli()
        self.assertEqual(record.step("replicates")["status"], COMPLETED)
        art = record.data["artifacts"]["feature_table_replicates"]
        # B ordered by injection: B3 (Rep1, failed), B2, B1 -> representative B2
        self.assertEqual(art["sample_columns"], ["A1", "B2", "D2", "E1", "QC1"])
        self.assertEqual(art["replicates"]["B2"], ["B2", "B1"])

        out = pd.read_csv(record.artifact_path("feature_table_replicates"), sep="\t")
        self.assertEqual(list(out.columns[:3]), ["id_number", "mz", "blank_flag"])
        np.testing.assert_allclose(
            out["A1"], self.ft[["A1", "A2", "A3"]].median(axis=1)
        )
        np.testing.assert_allclose(out["D2"], self.ft[["D2", "D3"]].median(axis=1))
        np.testing.assert_allclose(out["QC1"], self.ft["QC1"])

        path = record.artifact_path("replicate_correlations")
        self.assertEqual(pd.ExcelFile(path).sheet_names, ["replicate_correlations"])
        rep = pd.read_excel(path).set_index("biosample")
        self.assertEqual(list(rep.index), ["A", "B", "C", "D", "E"])
        for col in [
            "Rep1 vs Rep2",
            "Rep1 vs Rep3",
            "Rep2 vs Rep3",
            "min correlation",
            "max correlation",
            "average correlation",
            "flag",
        ]:
            self.assertIn(col, rep.columns)
        self.assertTrue(pd.isna(rep.at["A", "flag"]))
        self.assertEqual(rep.at["B", "flag"], "min < 0.7")
        self.assertEqual(rep.at["B", "Rep1"], "B3")
        self.assertEqual(rep.at["B", "dropped replicates"], "B3")
        self.assertEqual(rep.at["C", "flag"], "max < 0.7")
        self.assertEqual(rep.at["C", "status"], replicates.STATUS_EXCLUDED)
        self.assertTrue(pd.isna(rep.at["C", "Rep1 vs Rep3"]))

        summary = record.step("replicates")["summary"]
        self.assertEqual((summary["n_flagged"], summary["n_excluded"]), (3, 1))
        self.assertEqual(record.data["history"][-1]["status"], "completed")

    def test_no_drop_flagged(self):
        record = self.run_cli("--no_drop_flagged")
        art = record.data["artifacts"]["feature_table_replicates"]
        self.assertEqual(art["sample_columns"], ["A1", "B3", "C1", "D1", "E1", "QC1"])
        out = pd.read_csv(record.artifact_path("feature_table_replicates"), sep="\t")
        np.testing.assert_allclose(out["C1"], self.ft[["C1", "C2"]].median(axis=1))

    def test_cutoff_and_raw_intensities(self):
        record = self.run_cli("--pearson", "-1", "--no-log2")
        self.assertEqual(record.step("replicates")["summary"]["n_flagged"], 0)
        params = record.step("replicates")["params"]
        self.assertEqual((params["pearson"], params["no_log2"]), (-1, True))
        self.assertEqual(
            (params["rep_col"], params["table"]), ("biosample", "feature_table")
        )

    def test_metadata_option(self):
        other = self.meta.rename(columns={"biosample": "group_id"})
        path = self.root / "other_meta.csv"
        other.to_csv(path, index=False)
        cli.main(
            [
                "replicates",
                "-p",
                str(self.root),
                "--rep_col",
                "group_id",
                "--metadata",
                str(path),
            ]
        )
        record = ProjectRecord.load(str(self.root))
        self.assertEqual(record.step("replicates")["status"], COMPLETED)

    def test_missing_metadata(self):
        (self.root / "metadata.tsv").unlink()
        del self.record.data["artifacts"]["metadata"]
        self.record.save()
        with self.assertRaisesRegex(FileNotFoundError, "No sample metadata"):
            self.run_cli()
        record = ProjectRecord.load(str(self.root))
        self.assertEqual(record.step("replicates")["status"], FAILED)
        self.assertIn("No sample metadata", record.step("replicates")["message"])
        self.assertEqual(record.data["history"][-1]["status"], "failed")

    def test_missing_rep_col(self):
        with self.assertRaisesRegex(ValueError, "not in the metadata"):
            cli.main(["replicates", "-p", str(self.root), "--rep-col", "nope"])


if __name__ == "__main__":
    unittest.main()
