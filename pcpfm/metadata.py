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
import re
from collections.abc import Sequence
from typing import Any

import pandas as pd

MZML_EXT = ".mzml"

# sample_type words that mark non-study samples (blanks, QC, standards...)
NON_STUDY_TYPES = {
    "blank",
    "qc",
    "pool",
    "pooled",
    "qstd",
    "std",
    "standard",
    "nist",
    "dda",
    "solvent",
    "process",
}


def read_table(path: str, **kwargs: Any) -> pd.DataFrame:
    """Read a delimited table using its filename extension.

    Parameters
    ----------
    path : str
        File path to read or write.
    **kwargs : Any
        Additional keyword arguments forwarded to pandas.read_csv.

    Returns
    -------
    pd.DataFrame
        Table parsed by pandas.

    Notes
    -----
    Use tabs for .tsv, .txt and .tab; otherwise use commas.
    """
    sep = "\t" if path.lower().endswith((".tsv", ".txt", ".tab")) else ","
    return pd.read_csv(path, sep=sep, **kwargs)


def file_stem(path: str | os.PathLike[str]) -> str:
    """Extract the filename without a case-insensitive mzML suffix.

    Parameters
    ----------
    path : str | os.PathLike[str]
        File path to read or write.

    Returns
    -------
    str
        File stem with its original case preserved.
    """
    name = os.path.basename(str(path).strip())
    if name.lower().endswith(MZML_EXT):
        name = name[: -len(MZML_EXT)]
    return name


def load_metadata(path: str) -> pd.DataFrame:
    """Load and validate sample metadata.

    Parameters
    ----------
    path : str
        File path to read or write.

    Returns
    -------
    pd.DataFrame
        String-valued metadata with unique sample IDs and a sample_type column.

    Raises
    ------
    ValueError
        If the sample column is absent or contains empty or duplicate IDs.

    Notes
    -----
    Whitespace is stripped from column names and sample IDs. Missing sample
    types are empty strings.
    """
    df = read_table(path, dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    if "sample" not in df.columns:
        raise ValueError(
            f"Metadata {path} must have a 'sample' column; found {list(df.columns)}"
        )
    df["sample"] = df["sample"].str.strip()
    dup = df["sample"][df["sample"].duplicated()].tolist()
    if dup:
        raise ValueError(f"Duplicate sample IDs in metadata: {dup[:10]}")
    if (df["sample"] == "").any():
        raise ValueError("Empty sample ID in metadata.")
    if "sample_type" not in df.columns:
        df["sample_type"] = ""
    return df


def is_blank(sample_type: object) -> bool:
    """Identify sample types containing the word fragment blank.

    Parameters
    ----------
    sample_type : object
        Sample-type label to inspect after conversion to text.

    Returns
    -------
    bool
        Whether the label contains "blank", ignoring case.
    """
    return "blank" in str(sample_type).lower()


def sample_type_tokens(text: object) -> list[str]:
    """Split a sample type at separators, case changes and digit boundaries.

    Parameters
    ----------
    text : object
        Text value to parse.

    Returns
    -------
    list[str]
        Lowercase tokens, for example ["pooled", "qc"] for "pooledQC".
    """
    return [
        t.lower() for t in re.findall(r"[A-Z]+(?![a-z])|[A-Z]?[a-z]+|\d+", str(text))
    ]


def is_study_sample(sample_type: object) -> bool:
    """Identify study samples from their sample-type tokens.

    Parameters
    ----------
    sample_type : object
        Sample-type label to inspect after conversion to text.

    Returns
    -------
    bool
        False when a token, or its singular form, identifies a non-study sample.

    Notes
    -----
    Empty and unspecified labels are study samples. Labels such as QCpool lack a
    clear word boundary and may be classified as study samples; use an explicit
    sample selector for these labels.
    """
    for t in sample_type_tokens(sample_type):
        if t in NON_STUDY_TYPES or (t.endswith("s") and t[:-1] in NON_STUDY_TYPES):
            return False
    return True


def list_input_files(
    files_list: str | None = None, input_dirs: Sequence[str] | None = None
) -> list[str]:
    """Collect unique absolute mzML paths from a list and directories.

    Parameters
    ----------
    files_list : str | None, optional
        Optional text file containing one mzML path per line.
    input_dirs : Sequence[str] | None, optional
        Directories to scan, non-recursively, in the supplied order.

    Returns
    -------
    list[str]
        Paths in first-occurrence order: list entries followed by sorted
        directory entries.

    Notes
    -----
    Blank lines and comments beginning with # are ignored. Relative list entries
    resolve against the list file directory. Paths are normalized with abspath,
    without resolving symlinks.
    """
    paths = []
    if files_list:
        base = os.path.dirname(os.path.abspath(files_list))
        with open(files_list, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                if line and line.lower().endswith(MZML_EXT):
                    paths.append(
                        line if os.path.isabs(line) else os.path.join(base, line)
                    )
    for d in input_dirs or []:
        paths += sorted(
            p for p in glob.glob(os.path.join(d, "*")) if p.lower().endswith(MZML_EXT)
        )
    # normalize before de-duplicating: one file can be reached by an absolute path
    # in --files and by a relative --input directory
    return list(dict.fromkeys(os.path.abspath(p) for p in paths))


def resolve_inputs(
    metadata: pd.DataFrame | None,
    candidate_files: Sequence[str],
    metadata_path: str | None = None,
) -> tuple[list[dict[str, Any]], list[str]]:
    """Match metadata samples to mzML files.

    Parameters
    ----------
    metadata : pandas.DataFrame | None
        Sample metadata, or None to infer samples from candidate filenames.
    candidate_files : Sequence[str]
        Candidate absolute mzML paths.
    metadata_path : str | None, optional
        Metadata filename used to resolve relative file_path entries.

    Returns
    -------
    tuple[list[dict[str, Any]], list[str]]
        Input records and candidate files not referenced by metadata. Each
        record contains sample, sample_type, file, stem and blank; unresolved
        files are None.

    Raises
    ------
    ValueError
        If files share a stem or multiple metadata rows reference one file.

    Notes
    -----
    Without metadata, each candidate file becomes a sample named after its stem.
    """
    by_stem = {}
    for p in candidate_files:
        by_stem.setdefault(file_stem(p), []).append(p)
    clashes = {s: ps for s, ps in by_stem.items() if len(ps) > 1}
    if clashes:
        # Asari names sample columns by file stem, so stems must be unique
        raise ValueError(
            f"Different files share the same name, which Asari cannot distinguish: {clashes}"
        )

    if metadata is None:
        return [
            {"sample": s, "sample_type": "", "file": ps[0], "stem": s, "blank": False}
            for s, ps in by_stem.items()
        ], []

    base = (
        os.path.dirname(os.path.abspath(metadata_path))
        if metadata_path
        else os.getcwd()
    )
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
    # one file used by several rows; checked first as the more specific diagnosis
    by_path = {}
    for x in inputs:
        if x["file"]:
            by_path.setdefault(x["file"], []).append(x["sample"])
    shared = {p: s for p, s in by_path.items() if len(s) > 1}
    if shared:
        raise ValueError(
            f"Several metadata rows point to the same data file: {dict(list(shared.items())[:10])}"
        )
    # different files with one name; Asari names sample columns by file stem
    paths_by_stem = {}
    for p in by_path:
        paths_by_stem.setdefault(file_stem(p), []).append(p)
    clashes = {s: ps for s, ps in paths_by_stem.items() if len(ps) > 1}
    if clashes:
        raise ValueError(
            "Different files share the same name, which Asari cannot distinguish: "
            f"{dict(list(clashes.items())[:10])}"
        )
    unlisted = [p for p in candidate_files if p not in used]
    return inputs, unlisted
