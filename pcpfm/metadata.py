"""
Sample metadata and input file resolution.

Metadata is a CSV or TSV table with one row per sample (injection). Recognized columns:

    sample              required, unique sample ID; becomes the column name in feature tables
    original_file_name  data file name, with or without .mzML; used to match files
    file_path           optional full path to the mzML file (relative to the metadata file)
    sample_type         e.g. blank, qstd, pooled, sample; rows whose type contains "blank"
                        are processed in a separate Asari run
    group, batch, injection_order, dilution_factor
                        optional; used by QC when present

Other columns are carried along unchanged.
"""

import glob
import os

import pandas as pd

MZML_EXT = ".mzml"


def read_table(path, **kwargs):
    """Read CSV or TSV by extension (tab for .tsv/.txt)."""
    sep = "\t" if path.lower().endswith((".tsv", ".txt", ".tab")) else ","
    return pd.read_csv(path, sep=sep, **kwargs)


def file_stem(path):
    name = os.path.basename(str(path).strip())
    if name.lower().endswith(MZML_EXT):
        name = name[: -len(MZML_EXT)]
    return name


def load_metadata(path):
    df = read_table(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    if "sample" not in df.columns:
        raise ValueError(f"Metadata {path} must have a 'sample' column; found {list(df.columns)}")
    df["sample"] = df["sample"].str.strip()
    dup = df["sample"][df["sample"].duplicated()].tolist()
    if dup:
        raise ValueError(f"Duplicate sample IDs in metadata: {dup[:10]}")
    if (df["sample"] == "").any():
        raise ValueError("Empty sample ID in metadata.")
    if "sample_type" not in df.columns:
        df["sample_type"] = ""
    return df


def is_blank(sample_type):
    return "blank" in str(sample_type).lower()


def list_input_files(files_list=None, input_dirs=None):
    """Collect mzML paths from a --files list and/or --input directories."""
    paths = []
    if files_list:
        base = os.path.dirname(os.path.abspath(files_list))
        with open(files_list, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and line.lower().endswith(MZML_EXT):
                    paths.append(line if os.path.isabs(line) else os.path.join(base, line))
    for d in input_dirs or []:
        paths += sorted(
            p for p in glob.glob(os.path.join(d, "*")) if p.lower().endswith(MZML_EXT)
        )
    return [os.path.abspath(p) for p in dict.fromkeys(paths)]


def resolve_inputs(metadata, candidate_files, metadata_path=None):
    """
    Match samples to mzML files.

    Returns (inputs, unlisted_files):
        inputs: list of dicts {sample, sample_type, file (or None), stem, blank}
        unlisted_files: candidate files not referenced by the metadata

    Without metadata, every candidate file becomes a sample named by its file stem.
    """
    by_stem = {}
    for p in candidate_files:
        by_stem.setdefault(file_stem(p), []).append(p)
    clashes = {s: ps for s, ps in by_stem.items() if len(ps) > 1}
    if clashes:
        # Asari names sample columns by file stem, so stems must be unique
        raise ValueError(f"Different files share the same name, which Asari cannot distinguish: {clashes}")

    if metadata is None:
        return [
            {"sample": s, "sample_type": "", "file": ps[0], "stem": s, "blank": False}
            for s, ps in by_stem.items()
        ], []

    base = os.path.dirname(os.path.abspath(metadata_path)) if metadata_path else os.getcwd()
    inputs, used = [], set()
    for row in metadata.to_dict("records"):
        path = None
        if row.get("file_path"):
            p = row["file_path"]
            p = p if os.path.isabs(p) else os.path.join(base, p)
            path = os.path.abspath(p) if os.path.exists(p) else None
        else:
            for key in (row.get("original_file_name"), row["sample"]):
                if key and file_stem(key) in by_stem:
                    path = by_stem[file_stem(key)][0]
                    break
        if path:
            used.add(path)
        inputs.append(
            {
                "sample": row["sample"],
                "sample_type": row.get("sample_type", ""),
                "file": path,
                "stem": file_stem(path) if path else None,
                "blank": is_blank(row.get("sample_type", "")),
            }
        )
    stems = [x["stem"] for x in inputs if x["stem"]]
    if len(stems) != len(set(stems)):
        raise ValueError("Several metadata rows point to the same data file.")
    unlisted = [p for p in candidate_files if p not in used]
    return inputs, unlisted
