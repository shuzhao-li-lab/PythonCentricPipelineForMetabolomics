"""Tests for ``pcpfm.annotate.gc_references``. Run with `python -m unittest tests.test_annotate_gc_refs -v`."""

import os
import tempfile
import unittest

from pcpfm.annotate.gc_references import load_kovats


class TestKovats(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name

    def write(self, text):
        path = os.path.join(self.dir, "kovats.tsv")
        with open(path, "w") as f:
            f.write(text)
        return path

    def test_returns_rt_range_in_seconds(self):
        path = self.write("RI\tRT\n1000\t5.7\n1100\t6.3\n1200\t7.0\n")
        self.assertEqual(load_kovats(path), (5.7 * 60, 7.0 * 60))

    def test_rejects_unusable_files(self):
        bad = {
            "too few rows": "RI\tRT\n1000\t5.7\n",
            "one column": "RI\n1000\n1100\n",
            "non-numeric": "RI\tRT\n1000\t5.7\nC11\t6.3\n",
            "unsorted RT": "RI\tRT\n1000\t6.3\n1100\t5.7\n",
            "repeated RT": "RI\tRT\n1000\t5.7\n1100\t5.7\n",
            "decreasing RI": "RI\tRT\n1100\t5.7\n1000\t6.3\n",
            "infinite RT": "RI\tRT\n1000\t5.7\n1100\tinf\n",
            "no header": "1000\t5.7\n1100\t6.3\n1200\t7.0\n",
            "negative RT": "RI\tRT\n1000\t-2\n1100\t0\n",
            "columns swapped": "RT\tRI\n5.7\t1000\n6.3\t1100\n",
        }
        for case, text in bad.items():
            with self.subTest(case), self.assertRaises(ValueError):
                load_kovats(self.write(text))
