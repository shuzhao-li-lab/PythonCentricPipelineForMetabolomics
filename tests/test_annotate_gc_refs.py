"""Tests for ``pcpfm.annotate.gc_references``. Run with `python -m unittest tests.test_annotate_gc_refs -v`."""

import os
import tempfile
import unittest

from pcpfm.annotate.gc_references import check_db, load_kovats


def msp_entry(
    name, ri="1305", mass="116.095", peaks=("73.0396 999.0", "44.0131 412.75")
):
    """One MSP record; pass None to leave a field out."""
    lines = [f"Name: {name}"]
    if ri is not None:
        lines.append(f"RETENTIONTIME: {ri}")
    if mass is not None:
        lines.append(f"ExactMass: {mass}")
    if peaks is not None:
        lines += [f"Num Peaks: {len(peaks)}", *peaks]
    return "\n".join(lines) + "\n\n"


class TestCheckDb(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = tmp.name

    def write(self, *entries, name="lib.msp"):
        path = os.path.join(self.dir, name)
        with open(path, "w") as f:
            f.write("".join(entries))
        return path

    def test_counts_used_entries_and_warns_on_drops(self):
        path = self.write(msp_entry("A"), msp_entry("B"), msp_entry("no RI", ri=None))
        n, warnings = check_db(path)
        self.assertEqual(n, 2)
        self.assertEqual(len(warnings), 1)
        self.assertIn("1 entries", warnings[0])

    def test_entry_asari_cannot_read_fails(self):
        for case, entry in {
            "no ExactMass": msp_entry("A", mass=None),
            "text RI": msp_entry("A", ri="n/a"),
            "zero peaks": msp_entry("A", peaks=()),
            "nan RI": msp_entry("A", ri="nan"),
        }.items():
            with self.subTest(case), self.assertRaisesRegex(ValueError, "numeric"):
                check_db(self.write(entry))

    def test_no_usable_entry_fails(self):
        with self.assertRaisesRegex(ValueError, "no entry"):
            check_db(self.write(msp_entry("A", ri=None)))

    def test_other_format_fails(self):
        for name in ("lib.txt", "lib.Msp"):
            with self.subTest(name), self.assertRaisesRegex(ValueError, "MSP or JSON"):
                check_db(self.write(msp_entry("A"), name=name))


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
