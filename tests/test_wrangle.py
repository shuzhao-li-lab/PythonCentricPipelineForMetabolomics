"""Batch correction contract and numerical regression tests."""

import argparse
import contextlib
import io
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from pcpfm import cli, wrangle


class TestBatchCorrect(unittest.TestCase):
    def setUp(self):
        self.X = pd.DataFrame(
            {
                "a": [100.0, 200.0, 0.0],
                "b": [200.0, 400.0, 50.0],
                "c": [100.0, 200.0, 10.0],
                "d": [400.0, 800.0, 20.0],
            }
        )

    def test_batch_correct_reduces_batch_shift(self):
        feature = np.arange(20, dtype=float)
        baseline = 100 + 10 * feature
        shift = 25 + 3 * (feature % 5)
        spread1 = feature % 3 + 1
        spread2 = 1.5 * (feature % 4 + 1)
        X = pd.DataFrame(
            {
                "a": baseline - spread1,
                "b": baseline + spread1,
                "c": baseline + shift - spread2,
                "d": baseline + shift + spread2,
            },
            index=[f"F{i}" for i in range(20)],
        )
        original = X.copy()
        # Metadata order differs from table order; unrelated samples must be ignored.
        batches = pd.Series(
            {"d": "two", "unused": "three", "a": "one", "c": "two", "b": "one"}
        )
        corrected = wrangle.batch_correct(X, method="combat", batches=batches)
        before = (X[["c", "d"]].mean(axis=1) - X[["a", "b"]].mean(axis=1)).abs()
        after = (
            corrected[["c", "d"]].mean(axis=1) - corrected[["a", "b"]].mean(axis=1)
        ).abs()
        self.assertTrue((after < before * 0.1).all())
        self.assertTrue(np.isfinite(corrected.to_numpy()).all())
        np.testing.assert_allclose(
            corrected.mean(axis=1), X.mean(axis=1), atol=0.25, rtol=0
        )
        pd.testing.assert_frame_equal(X, original)
        pd.testing.assert_index_equal(corrected.index, X.index)
        pd.testing.assert_index_equal(corrected.columns, X.columns)

        args = argparse.Namespace(
            batch_field="batch",
            normalize="none",
            batch_correct="combat",
            impute="half-min",
            no_log2=False,
        )
        prepared, log = wrangle.prepare_analysis_table(
            X, batches.to_frame("batch"), args
        )
        pd.testing.assert_frame_equal(prepared, wrangle.log2(corrected))
        self.assertIn("batch correction: combat", log)

    def test_batch_correct_validation(self):
        self.assertIs(
            wrangle.batch_correct(self.X, method="none"),
            self.X,
        )
        for batches in (
            None,
            pd.Series({"a": 1, "b": 1, "c": 1, "d": 1}),
            pd.Series({"a": 1, "b": 1, "c": 2}),
            pd.Series({"a": "one", "b": "one", "c": "two", "d": " "}),
        ):
            with self.subTest(batches=batches), self.assertRaises(ValueError):
                wrangle.batch_correct(
                    self.X,
                    method="combat",
                    batches=batches,
                )
        with self.assertRaisesRegex(ValueError, "Unknown batch correction"):
            wrangle.batch_correct(self.X, method="unknown")
        batches = pd.Series({"a": 1, "b": 1, "c": 2, "d": 2})
        for value in (np.nan, np.inf):
            X = self.X.copy()
            X.at[0, "a"] = value
            with (
                self.subTest(value=value),
                self.assertRaises(ValueError),
            ):
                wrangle.batch_correct(X, method="combat", batches=batches)

    def test_batch_correct_nonfinite_output(self):
        batches = pd.Series({"a": 1, "b": 1, "c": 2, "d": 2})
        result = self.X.copy()
        result["a"] = [-2.0, 0.0, 0.5]
        with patch("pcpfm.wrangle.pycombat_norm", return_value=result):
            corrected = wrangle.batch_correct(
                self.X,
                method="combat",
                batches=batches,
            )
        expected = result.copy()
        expected["a"] = [-2.0, 0.0, 0.5]
        pd.testing.assert_frame_equal(corrected, expected)
        for value in (np.nan, np.inf, -np.inf):
            result.at[0, "a"] = value
            with self.subTest(value=value), patch(
                "pcpfm.wrangle.pycombat_norm", return_value=result
            ):
                with self.assertRaisesRegex(ValueError, "non-finite"):
                    wrangle.batch_correct(self.X, method="combat", batches=batches)

    def test_cli_choices(self):
        parser = cli.build_parser()
        default = parser.parse_args(["qc", "-p", "project"])
        self.assertEqual(default.batch_correct, "none")
        for method in ("none", "combat"):
            args = parser.parse_args(["qc", "-p", "project", "--batch-correct", method])
            self.assertEqual(args.batch_correct, method)
        with patch("pcpfm.wrangle.pycombat_norm") as combat:
            with (
                contextlib.redirect_stderr(io.StringIO()),
                self.assertRaises(SystemExit) as error,
            ):
                cli.main(["qc", "-p", "project", "--batch-correct", "invalid"])
            self.assertEqual(error.exception.code, 2)
            combat.assert_not_called()

    def test_bad_sample_values_and_duplicate_batch_labels(self):
        batches = pd.Series({"a": 1, "b": 1, "c": 2, "d": 2})
        bad = self.X.copy()
        bad["a"] = ["bad", "input", "values"]
        # InMoose rejects nonnumeric input; the wrapper propagates its TypeError.
        with self.assertRaises(TypeError):
            wrangle.batch_correct(bad, method="combat", batches=batches)
        with self.assertRaisesRegex(ValueError, "unique sample names"):
            wrangle.batch_correct(
                self.X,
                method="combat",
                batches=pd.Series([1, 2], index=["a", "a"]),
            )

    def test_empty_table_is_rejected_before_combat(self):
        for table in (self.X.iloc[:0], self.X.iloc[:, :0]):
            with self.subTest(shape=table.shape), patch(
                "pcpfm.wrangle.pycombat_norm"
            ) as combat:
                with self.assertRaisesRegex(ValueError, "input table is empty"):
                    wrangle.batch_correct(table, method="combat")
                combat.assert_not_called()

    def test_known_input_matches_inmoose_reference(self):
        # Frozen output from a direct InMoose 0.9.1 pycombat_norm call.
        table = pd.DataFrame(
            {
                "a": [99.0, 108.0, 117.0, 129.0],
                "b": [101.0, 112.0, 123.0, 131.0],
                "c": [123.5, 135.0, 146.5, 157.0],
                "d": [126.5, 141.0, 155.5, 169.0],
            },
            index=["F1", "F2", "F3", "F4"],
        )
        columns = list(table.columns)
        original = table.copy()
        batches = pd.Series({"d": 2, "a": 1, "c": 2, "b": 1})
        expected = pd.DataFrame(
            [
                [
                    110.78977498293092,
                    113.6069156936462,
                    111.7699330288871,
                    114.19777109666532,
                ],
                [
                    121.21502540919742,
                    126.87645344972046,
                    121.49434935732272,
                    126.3584552423661,
                ],
                [
                    131.64974848497167,
                    140.13445205067168,
                    131.22252286820068,
                    138.5163569096646,
                ],
                [
                    145.37467215738874,
                    148.5233264220144,
                    140.85923168143776,
                    150.54039117447357,
                ],
            ],
            index=table.index,
            columns=columns,
        )
        actual = wrangle.batch_correct(table, method="combat", batches=batches)
        pd.testing.assert_frame_equal(actual[columns], expected, rtol=1e-7, atol=1e-7)
        pd.testing.assert_frame_equal(table, original)
