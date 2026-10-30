"""
Evidence table and manifest for ``pcpfm annotate lc``.

The evidence table is long: one row per accepted candidate, never one row per
feature, so it can hold several candidates for one feature without touching the
quantitative table. Columns an evidence type does not measure stay empty.
The manifest is the provenance record of one annotation run.
"""

import hashlib
import json
import os
from typing import Any

import pandas as pd

SCHEME = "pcpfm_v1"

EVIDENCE_COLUMNS = [
    "feature_id",
    "empcpd_id",
    "spectrum_id",
    "scheme",
    "level",
    "evidence_type",
    "reference_id",
    "reference_name",
    "reference_source",
    "neutral_formula",
    "ion_relation",
    "mz_error_ppm",
    "rt_error_s",
    "score",
    "n_matched_peaks",
]


def empty_evidence() -> pd.DataFrame:
    """Return an evidence table with all columns and no rows.

    Returns
    -------
    pd.DataFrame
        Zero-row table with ``EVIDENCE_COLUMNS``.
    """
    return pd.DataFrame(columns=EVIDENCE_COLUMNS)


def write_evidence(rows: pd.DataFrame, path: str) -> int:
    """Write evidence rows as TSV in the fixed column order.

    Parameters
    ----------
    rows : pd.DataFrame
        Evidence rows; missing columns are written empty.
    path : str
        Output TSV path.

    Returns
    -------
    int
        Number of rows written; zero writes a header-only file.

    Raises
    ------
    ValueError
        If ``rows`` has columns outside ``EVIDENCE_COLUMNS``.
    """
    unknown = [c for c in rows.columns if c not in EVIDENCE_COLUMNS]
    if unknown:
        raise ValueError(f"Unknown evidence columns: {unknown}")
    rows.reindex(columns=EVIDENCE_COLUMNS).to_csv(path, sep="\t", index=False)
    return len(rows)


def file_sha256(path: str) -> str:
    """Return the SHA-256 hex digest of a file's content.

    Parameters
    ----------
    path : str
        File to hash.

    Returns
    -------
    str
        Hex digest.
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_manifest(manifest: dict[str, Any], path: str) -> None:
    """Write the run manifest as indented JSON, replacing any old file atomically.

    Parameters
    ----------
    manifest : dict[str, Any]
        JSON-serializable provenance record.
    path : str
        Output JSON path.
    """
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, default=str)
    os.replace(tmp, path)
