"""
Helper module to resolve cyclic imports
Provided feature-table loading functions and sample metadata shared by pipeline steps.
"""

import os

import pandas as pd

FEATURE_COLUMNS = [
    "id_number",
    "mz",
    "rtime",
    "rtime_left_base",
    "rtime_right_base",
    "parent_masstrack_id",
    "peak_area",
    "cSelectivity",
    "goodness_fitting",
    "snr",
    "detection_counts",
]


BLANK_COLUMNS = [
    "blank_match", "blank_mz", "blank_rtime", "blank_ppm", "blank_drt",
    "blank_mean", "sample_mean", "blank_ratio", "blank_flag",
]


def _validate_samples(ft, columns):
    if not columns or len(set(columns)) != len(columns):
        raise ValueError("A feature table needs a non-empty list of unique sample columns.")
    missing = [c for c in columns if c not in ft.columns]
    if missing:
        raise ValueError(f"Recorded sample columns are missing from the table: {missing}")


def load_asari_table(path, stem_to_sample):
    """
    Read Asari output, identifying samples only by resolved input file stems.

    Sample columns are renamed to sample IDs and kept in table order. Any other column
    Asari adds in a future version is kept as a feature annotation and reported once,
    never guessed to be a sample.
    """
    ft = pd.read_csv(path, sep="\t", dtype={"id_number": str})
    columns = [c for c in ft.columns if c in stem_to_sample]
    unknown = [c for c in ft.columns if c not in stem_to_sample and c not in FEATURE_COLUMNS]
    if unknown:
        print(f"[pcpfm] Note: Asari table has columns pcpfm does not know: {unknown}; "
              "kept as feature annotations.")
    _validate_samples(ft, columns)
    ft = ft.rename(columns={c: stem_to_sample[c] for c in columns})
    if not ft.columns.is_unique:
        raise ValueError("Renaming Asari samples would create duplicate table columns.")
    return ft, [stem_to_sample[c] for c in columns]


def load_table(record, name):
    """Read a table using its artifact's explicit sample columns.

    Older artifacts use known project input sample names present in the table.
    Unknown columns are retained as metadata, never inferred to be samples.
    """
    ft = pd.read_csv(record.artifact_path(name), sep="\t", dtype={"id_number": str})
    columns = record.data["artifacts"][name].get("sample_columns")
    if columns is None:
        # artifacts from an earlier dev version have no sample list
        print(f"[pcpfm] Warning: artifact '{name}' has no sample list; inferring sample columns "
              "from the project inputs.")
        samples = {entry["sample"] for entry in record.data["inputs"]}
        columns = [c for c in ft.columns if c in samples]
        if not columns:
            raise ValueError(
                f"Artifact '{name}' has no sample_columns or matching project inputs; "
                "rerun process to register its sample columns."
            )
    if not isinstance(columns, list):
        raise ValueError(f"Artifact '{name}' sample_columns must be a list.")
    _validate_samples(ft, columns)
    return ft, list(columns)


def parse_selector(text):
    """'sample_type=pooled' -> ('sample_type', {'pooled'}); values separated by commas."""
    if not text:
        return None
    key, _, values = text.partition("=")
    return key.strip(), {v.strip().lower() for v in values.split(",") if v.strip()}


def sample_info(record, sample_cols):
    """Metadata rows for the table's samples, in table order, with injection_order filled."""
    path = record.artifact_path("metadata", must_exist=False)
    if path and os.path.exists(path):
        meta = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        meta["_order"] = range(len(meta))
        meta = meta.set_index("sample").reindex(sample_cols)
    else:
        meta = pd.DataFrame(index=pd.Index(sample_cols, name="sample"))
        meta["_order"] = range(len(meta))
    meta["sample_type"] = meta.get(
        "sample_type", pd.Series("", index=meta.index)
    ).fillna("")
    meta.loc[meta["sample_type"] == "", "sample_type"] = "unspecified"
    if "injection_order" in meta.columns:
        meta["injection_order"] = pd.to_numeric(
            meta["injection_order"], errors="coerce"
        )
    else:
        meta["injection_order"] = meta["_order"]
    return meta.drop(columns="_order")


def select(meta, selector):
    if selector is None:
        return []
    key, values = selector
    if key not in meta.columns:
        return []
    return [s for s in meta.index if str(meta.at[s, key]).lower() in values]
