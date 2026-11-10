#!/usr/bin/env python3
"""Summarize mzML run headers into a TSV, one row per file.

Only the header and the first spectra are read, so large files are fast.
Standard library only, so it runs without the pcpfm environment.

Usage: python mzml_summary.py DIR_OR_FILE [DIR_OR_FILE ...] [-o summary.tsv]
Directories are searched recursively.
"""

import argparse
import csv
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

N_SPECTRA = 50  # spectra inspected for polarity and MS levels
COLUMNS = [
    "file",
    "stem",
    "start_time",
    "acq_rank",  # order of start_time within the same ionization and polarity
    "ionization",
    "polarity",
    "ms_levels",
    "n_spectra",
]


def summarize(path: Path) -> dict:
    """Read acquisition facts from one mzML file.

    Parameters
    ----------
    path : Path
        mzML file to read.

    Returns
    -------
    dict
        Values for COLUMNS except acq_rank; unreadable fields stay empty.
    """
    name = path.name
    stem = name[:-5] if name.lower().endswith(".mzml") else path.stem
    row = dict.fromkeys(COLUMNS, "") | {"file": str(path), "stem": stem}
    ionization, polarity, levels = [], set(), set()
    in_source, seen = False, 0
    try:
        for event, el in ET.iterparse(path, events=("start", "end")):
            tag = el.tag.rsplit("}", 1)[-1]
            if event == "start":
                if tag == "source":
                    in_source = True
                elif tag == "run":
                    row["start_time"] = el.get("startTimeStamp", "")
                elif tag == "spectrumList":
                    row["n_spectra"] = el.get("count", "")
                continue
            if tag == "source":
                in_source = False
            elif tag == "cvParam":
                cv = el.get("name", "")
                if in_source and cv.endswith("ionization") and cv not in ionization:
                    ionization.append(cv)
                elif cv in ("positive scan", "negative scan"):
                    polarity.add(cv.split()[0][:3])
                elif cv == "ms level":
                    levels.add(int(el.get("value", 0)))
            elif tag == "spectrum":
                seen += 1
                el.clear()
                if seen >= N_SPECTRA:
                    break
    except (ET.ParseError, OSError, ValueError) as e:
        print(f"Warning: could not fully read {path}: {e}", file=sys.stderr)
    row["ionization"] = ";".join(ionization)
    row["polarity"] = "+".join(sorted(polarity))
    row["ms_levels"] = ",".join(str(x) for x in sorted(levels))
    return row


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("paths", nargs="+", help="mzML files or directories")
    ap.add_argument("-o", "--output", help="output TSV (default: stdout)")
    args = ap.parse_args(argv)

    files = []
    for p in map(Path, args.paths):
        if p.is_dir():
            files += sorted(
                f for f in p.rglob("*")
                if f.suffix.lower() == ".mzml" and not f.name.startswith(".")
            )
        else:
            files.append(p)
    rows = [summarize(f) for f in files]

    # acquisition order within each mode; files without a timestamp get no rank
    modes = {}
    for r in sorted((r for r in rows if r["start_time"]), key=lambda r: r["start_time"]):
        mode = modes.setdefault((r["ionization"], r["polarity"]), [])
        mode.append(r)
        r["acq_rank"] = len(mode)

    stems = [r["stem"] for r in rows]
    clashes = sorted({s for s in stems if stems.count(s) > 1})
    if clashes:
        print(f"Warning: files in different folders share a name: {clashes}", file=sys.stderr)

    out = open(args.output, "w", newline="", encoding="utf-8") if args.output else sys.stdout
    writer = csv.DictWriter(out, COLUMNS, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    if args.output:
        out.close()


if __name__ == "__main__":
    main()
