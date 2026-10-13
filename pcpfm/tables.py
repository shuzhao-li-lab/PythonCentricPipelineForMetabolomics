"""Load feature tables and align sample metadata for shared pipeline steps."""

import os

import pandas as pd

from .project import ProjectRecord

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
    "blank_match",
    "blank_mz",
    "blank_rtime",
    "blank_ppm",
    "blank_drt",
    "blank_mean",
    "sample_mean",
    "blank_ratio",
    "blank_flag",
]


def _validate_samples(ft: pd.DataFrame, columns: list[str]) -> None:
    """Validate an explicit list of feature-table sample columns.

    Parameters
    ----------
    ft : pandas.DataFrame
        Feature table containing sample intensities and available feature
        annotations.
    columns : list[str]
        Sample column names to validate against the feature table.

    Raises
    ------
    ValueError
        If the list is empty, contains duplicates, or names missing columns.
    """
    if not columns or len(set(columns)) != len(columns):
        raise ValueError(
            "A feature table needs a non-empty list of unique sample columns."
        )
    missing = [c for c in columns if c not in ft.columns]
    if missing:
        raise ValueError(
            f"Recorded sample columns are missing from the table: {missing}"
        )


def load_asari_table(
    path: str | os.PathLike[str], stem_to_sample: dict[str, str]
) -> tuple[pd.DataFrame, list[str]]:
    """Load an Asari table and rename recognized sample columns.

    Parameters
    ----------
    path : str | os.PathLike[str]
        File path to read or write.
    stem_to_sample : dict[str, str]
        Mapping from resolved mzML file stems to sample IDs.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        Feature table and sample IDs in table-column order.

    Raises
    ------
    ValueError
        If no sample columns match or renaming creates duplicate columns.

    Notes
    -----
    Only resolved input stems identify samples. Unknown columns remain feature
    annotations and are reported. Feature metadata columns are optional.
    """
    ft = pd.read_csv(path, sep="\t", dtype={"id_number": str})
    columns = [c for c in ft.columns if c in stem_to_sample]
    unknown = [
        c for c in ft.columns if c not in stem_to_sample and c not in FEATURE_COLUMNS
    ]
    if unknown:
        print(
            f"[pcpfm] Note: Asari table has columns pcpfm does not know: {unknown}; "
            "kept as feature annotations."
        )
    _validate_samples(ft, columns)
    ft = ft.rename(columns={c: stem_to_sample[c] for c in columns})
    if not ft.columns.is_unique:
        raise ValueError("Renaming Asari samples would create duplicate table columns.")
    return ft, [stem_to_sample[c] for c in columns]


def load_table(record: ProjectRecord, name: str) -> tuple[pd.DataFrame, list[str]]:
    """Load a feature-table artifact using its recorded sample columns.

    Parameters
    ----------
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.
    name : str
        Feature-table artifact name.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        Feature table and an independent list of sample columns in recorded
        order.

    Raises
    ------
    TypeError
        If recorded sample_columns is not a list.
    ValueError
        If sample columns are empty, duplicated, or absent from the table.

    Notes
    -----
    Legacy artifacts without sample_columns use known project input names.
    Unknown columns remain annotations; metadata columns are not required.
    """
    ft = pd.read_csv(record.artifact_path(name), sep="\t", dtype={"id_number": str})
    columns = record.data["artifacts"][name].get("sample_columns")
    if columns is None:
        # artifacts from an earlier dev version have no sample list
        print(
            f"[pcpfm] Warning: artifact '{name}' has no sample list; inferring sample columns "
            "from the project inputs."
        )
        samples = {entry["sample"] for entry in record.data["inputs"]}
        columns = [c for c in ft.columns if c in samples]
        if not columns:
            raise ValueError(
                f"Artifact '{name}' has no sample_columns or matching project inputs; "
                "rerun process to register its sample columns."
            )
    if not isinstance(columns, list):
        raise TypeError(f"Artifact '{name}' sample_columns must be a list.")
    _validate_samples(ft, columns)
    return ft, list(columns)


def parse_selector(text: str | None) -> tuple[str, set[str]] | None:
    """Parse a column=value selector with comma-separated values.

    Parameters
    ----------
    text : str | None
        Selector such as "sample_type=pooled,qc", or None.

    Returns
    -------
    tuple[str, set[str]] | None
        Column name and stripped lowercase values, or None for empty input.
    """
    if not text:
        return None
    key, _, values = text.partition("=")
    return key.strip(), {v.strip().lower() for v in values.split(",") if v.strip()}


def sample_info(record: ProjectRecord, sample_cols: list[str]) -> pd.DataFrame:
    """Align sample metadata to feature-table sample columns.

    Parameters
    ----------
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.
    sample_cols : list[str]
        Sample intensity column names in the feature table, in the desired
        order.

    Returns
    -------
    pd.DataFrame
        Metadata indexed by sample ID in table order, with sample_type and
        injection_order.

    Notes
    -----
    Missing sample types become unspecified. If injection_order is absent,
    metadata row order is used, or table order when no metadata artifact exists.
    Supplied nonnumeric injection orders become NaN.
    """
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


def select(meta: pd.DataFrame, selector: tuple[str, set[str]] | None) -> list[str]:
    """Select sample IDs by a case-insensitive metadata value match.

    Parameters
    ----------
    meta : pandas.DataFrame
        Sample metadata indexed by sample ID.
    selector : tuple[str, set[str]] | None
        Column and lowercase allowed values, as returned by parse_selector.

    Returns
    -------
    list[str]
        Matching sample IDs in metadata order; empty if the selector or column
        is absent.
    """
    if selector is None:
        return []
    key, values = selector
    if key not in meta.columns:
        return []
    return [s for s in meta.index if str(meta.at[s, key]).lower() in values]
