"""Tests for ``pcpfm annotate gc``. Run with `python -m unittest tests.test_annotate_gc -v`."""

import contextlib
import io
import unittest

from pcpfm import cli

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
