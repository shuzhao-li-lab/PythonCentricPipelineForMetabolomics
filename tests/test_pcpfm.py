"""Unit tests for pcpfm v2. Run with `python -m unittest discover tests` or `pytest tests`."""

import argparse
import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from pcpfm import analyze, asari_runner, blanks, metadata, qc, wrangle
from pcpfm.project import COMPLETED, ProjectRecord


def feature_table(mz, rt, intensities, prefix="F"):
    ft = pd.DataFrame({"id_number": [f"{prefix}{i}" for i in range(len(mz))], "mz": mz, "rtime": rt})
    for name, values in intensities.items():
        ft[name] = values
    return ft


class TestBlanks(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(0)
        n = 200
        self.mz = np.sort(rng.uniform(80, 900, n))
        self.rt = rng.uniform(10, 300, n)
        self.samples = feature_table(self.mz, self.rt, {"s1": np.full(n, 1000.0), "s2": np.full(n, 3000.0)})
        # blank run: same compounds, m/z within 1 ppm, rtime shifted by +4 s; first 50 are strong in blanks
        blank_int = np.where(np.arange(n) < 50, 1500.0, 100.0)
        self.blank = feature_table(self.mz * (1 + 1e-6), self.rt + 4, {"b1": blank_int, "b2": blank_int}, prefix="B")

    def test_rt_offset(self):
        offset, pairs = blanks.estimate_rt_offset(
            self.mz, self.rt, self.blank["mz"].to_numpy(), self.blank["rtime"].to_numpy())
        self.assertAlmostEqual(offset, 4, places=3)
        self.assertGreater(pairs, 100)

    def test_match_and_flag(self):
        ft, summary = blanks.blank_filter_separate_runs(self.samples, ["s1", "s2"], self.blank, ["b1", "b2"])
        self.assertEqual(summary["n_matched"], 200)
        self.assertTrue((ft["blank_match"].str[1:] == ft["id_number"].str[1:]).all())
        # sample mean 2000; blank 1500 * 3 >= 2000 -> flagged; 100 * 3 < 2000 -> kept
        self.assertEqual(summary["n_flagged"], 50)
        self.assertTrue(ft["blank_flag"].iloc[:50].all())

    def test_no_match_outside_tolerance(self):
        far = self.blank.copy()
        far["rtime"] += 60
        m = blanks.match_features(self.samples, far, ppm=5, rt_tol=10)
        self.assertTrue((m["blank_match"] == "").all())

    def test_closest_candidate(self):
        s = feature_table([200.0], [100.0], {"s1": [10.0]})
        b = feature_table([200.0005, 200.0001], [103.0, 101.0], {"b1": [1.0, 2.0]}, prefix="B")
        m = blanks.match_features(s, b, ppm=5, rt_tol=10)
        self.assertEqual(m["blank_match"].iloc[0], "B1")

    def test_same_run(self):
        ft = feature_table([100.0, 200.0], [10.0, 20.0], {"s1": [100.0, 1000.0], "bl": [50.0, 0.0]})
        out, summary = blanks.blank_filter_same_run(ft, ["s1"], ["bl"])
        self.assertEqual(list(out["blank_flag"]), [True, False])


class TestMetadata(unittest.TestCase):
    def test_resolve_inputs(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ("a_001.mzML", "b_002.mzML", "extra.mzML"):
                open(os.path.join(d, name), "w").close()
            meta_path = os.path.join(d, "meta.csv")
            pd.DataFrame({
                "sample": ["A", "B", "C"],
                "original_file_name": ["a_001", "b_002.mzML", "missing"],
                "sample_type": ["sample", "Blank", "sample"],
            }).to_csv(meta_path, index=False)
            meta = metadata.load_metadata(meta_path)
            files = metadata.list_input_files(input_dirs=[d])
            inputs, unlisted = metadata.resolve_inputs(meta, files, meta_path)
            by = {x["sample"]: x for x in inputs}
            self.assertEqual(by["A"]["stem"], "a_001")
            self.assertTrue(by["B"]["blank"])
            self.assertIsNone(by["C"]["file"])
            self.assertEqual([os.path.basename(p) for p in unlisted], ["extra.mzML"])

    def test_files_list_and_no_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "x.mzML"), "w").close()
            lst = os.path.join(d, "files.txt")
            with open(lst, "w") as f:
                f.write("x.mzML\n# comment\nnot_mzml.raw\n")
            files = metadata.list_input_files(files_list=lst)
            self.assertEqual(len(files), 1)
            self.assertTrue(os.path.samefile(files[0], os.path.join(d, "x.mzML")))
            inputs, _ = metadata.resolve_inputs(None, files)
            self.assertEqual(inputs[0]["sample"], "x")

    def test_duplicate_samples(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "m.tsv")
            pd.DataFrame({"sample": ["A", "A"]}).to_csv(p, sep="\t", index=False)
            with self.assertRaises(ValueError):
                metadata.load_metadata(p)


class TestAsariLog(unittest.TestCase):
    def test_parse_alignment(self):
        log = (
            "Building composite mass tracks and calibrating retention time ...\n\n"
            "    Blank_003\n\t\tgood_reference_landmark_peaks:  86\n"
            "    Blank_001\n\tgood_landmark_peaks:  50\n"
            "    Bad_sample\n\tgood_landmark_peaks:  3\n"
            "    ~warning~ Failure in retention time alignment (3); Bad_sample.\n"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".log", delete=False) as f:
            f.write(log)
        try:
            a = asari_runner.parse_alignment(f.name)
        finally:
            os.remove(f.name)
        self.assertTrue(a["Blank_003"]["is_reference"])
        self.assertEqual(a["Blank_001"]["landmark_peaks"], 50)
        self.assertFalse(a["Bad_sample"]["aligned"])


class TestQC(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1)
        n_feat, samples = 300, [f"S{i}" for i in range(12)]
        base = rng.lognormal(10, 1.5, n_feat)
        data = {s: base * rng.lognormal(0, 0.1, n_feat) for s in samples}
        data["S11"] = base * rng.lognormal(0, 1.5, n_feat)       # outlier sample
        self.ft = feature_table(np.linspace(100, 800, n_feat), rng.uniform(0, 300, n_feat), data)
        self.ft["blank_flag"] = np.arange(n_feat) < 10
        self.samples = samples
        self.meta = pd.DataFrame(
            {"sample_type": ["pooled"] * 4 + ["sample"] * 8,
             "injection_order": range(12),
             "dilution_factor": [""] * 4 + ["1", "1", "2", "2", "4", "4", "8", "8"]},
            index=pd.Index(samples, name="sample"))

    def test_metrics_and_filter(self):
        X = self.ft[self.samples]
        pca_result = qc.pca(X)
        smet, corr = qc.sample_metrics(X, self.meta, {}, pca_result)
        self.assertIn("low_correlation", smet.at["S11", "flags"])
        pooled = qc.select(self.meta, qc.parse_selector("sample_type=pooled"))
        self.assertEqual(pooled, ["S0", "S1", "S2", "S3"])
        fmet = qc.feature_metrics(self.ft, X, self.meta, pooled, [s for s in self.samples if s not in pooled])
        self.assertLess(fmet["qc_cv"].median(), 20)
        self.assertIn("dilution_r", fmet.columns)
        args = argparse.Namespace(keep_blank_features=False, min_detection=0.0, max_qc_cv=0,
                                  drop_flagged_samples=True)
        out, cols, log = qc.filter_table(self.ft, self.samples, fmet, smet, args)
        self.assertEqual(len(out), 290)
        self.assertNotIn("S11", cols)


class TestWrangle(unittest.TestCase):
    def setUp(self):
        self.X = pd.DataFrame({"a": [100.0, 200.0, 0.0], "b": [200.0, 400.0, 50.0],
                               "c": [100.0, 200.0, 10.0], "d": [400.0, 800.0, 20.0]})

    def test_tic(self):
        # TIC over features present in > 90% of samples (first two rows); equalized to the median
        N = wrangle.normalize(self.X, "tic")
        tics = N.iloc[:2].sum()
        self.assertTrue(np.allclose(tics, tics.iloc[0]))

    def test_tic_by_batch(self):
        # features used for TIC are chosen within each batch; make row 3 absent in both batches
        X = self.X.copy()
        X.at[2, "c"] = 0.0
        batches = pd.Series({"a": "1", "b": "1", "c": "2", "d": "2"})
        N = wrangle.normalize(X, "tic", batches)
        self.assertTrue(np.allclose(N.iloc[:2].sum(), N.iloc[:2].sum().iloc[0]))

    def test_impute_and_log(self):
        I = wrangle.impute(self.X)
        self.assertEqual(I.at[2, "a"], 5.0)
        self.assertTrue(np.allclose(wrangle.log2(pd.DataFrame({"x": [0.0, 1.0]}))["x"], [0, 1]))

    def test_batch_correct_placeholder(self):
        self.assertIs(wrangle.batch_correct(self.X, "none"), self.X)
        with self.assertRaises(NotImplementedError):
            wrangle.batch_correct(self.X, "combat", pd.Series({"a": 1, "b": 1, "c": 2, "d": 2}))


class TestAnalyze(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(2)
        self.samples = [f"S{i}" for i in range(12)]
        self.Y = pd.DataFrame(rng.normal(20, 1, (200, 12)), columns=self.samples)
        self.group = pd.Series(["ctrl"] * 6 + ["case"] * 6, index=self.samples)
        self.age = pd.Series(np.linspace(20, 70, 12), index=self.samples)
        self.Y.iloc[0, 6:] += 5                                 # differential feature
        self.Y.iloc[1] = 10 + 0.2 * self.age.to_numpy()          # age-dependent feature

    def test_classify(self):
        meta = pd.DataFrame({"group": self.group, "age": self.age.astype(str),
                             "id": self.samples, "sample_type": "sample"})
        v = analyze.select_variables(meta, self.samples)
        self.assertEqual(v["group"][0], "categorical")
        self.assertEqual(v["age"][0], "continuous")
        self.assertNotIn("id", v)

    def test_study_samples(self):
        meta = pd.DataFrame({"sample_type": ["sample", "blank", "pooled", "qstd", "Unknown"]},
                            index=list("abcde"))
        self.assertEqual(analyze.study_samples(meta), ["a", "e"])

    def test_two_groups(self):
        res, desc, _ = analyze.test_two_groups(self.Y, self.group, reference="ctrl")
        self.assertEqual(desc, "case vs ctrl")
        self.assertAlmostEqual(res["log2FC"].iloc[0], 5, delta=1)
        self.assertEqual(res["p"].idxmin(), 0)

    def test_anova(self):
        g = pd.Series(["a"] * 4 + ["b"] * 4 + ["c"] * 4, index=self.samples)
        res, _, _ = analyze.test_anova(self.Y, g)
        self.assertTrue(((res["eta_squared"] >= 0) & (res["eta_squared"] <= 1)).all())

    def test_regression(self):
        res, _, _ = analyze.test_regression(self.Y, self.age)
        self.assertAlmostEqual(res["slope"].iloc[1], 0.2, places=6)
        self.assertAlmostEqual(res["r"].iloc[1], 1.0, places=6)

    def test_bh(self):
        q = analyze.bh_fdr([0.01, 0.04, 0.03, np.nan])
        self.assertTrue(np.allclose(q[:3], [0.03, 0.04, 0.04]))
        self.assertTrue(np.isnan(q[3]))

    def test_natural_sort(self):
        self.assertEqual(sorted(["V(16:1)", "V(4:1)", "V(1:1)"], key=analyze.natural_key),
                         ["V(1:1)", "V(4:1)", "V(16:1)"])


class TestProjectRecord(unittest.TestCase):
    def test_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            r = ProjectRecord.load_or_create(d, "demo")
            r.start_step("process")
            path = os.path.join(r.subdir("tables"), "t.tsv")
            open(path, "w").close()
            r.add_artifact("feature_table", path, "test", "process")
            r.finish_step("process")
            r2 = ProjectRecord.load(d)
            self.assertEqual(r2.data["project"]["name"], "demo")
            self.assertEqual(r2.data["artifacts"]["feature_table"]["path"], os.path.join("tables", "t.tsv"))
            self.assertTrue(r2.has_artifact("feature_table"))
            self.assertEqual(r2.step("process")["status"], COMPLETED)


if __name__ == "__main__":
    unittest.main()
