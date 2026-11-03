"""
Evidence table and manifest for ``pcpfm annotate lc``.

The evidence table is long: one row per accepted candidate, never one row per
feature, so it can hold several candidates for one feature without touching the
quantitative table. Columns an evidence type does not measure stay empty.
The manifest is the provenance record of one annotation run. Runs are written
through :func:`staged_dir`, so a failed run never touches the previous outputs.
"""

import contextlib
import hashlib
import json
import os
import shutil
import tempfile
from collections.abc import Iterator
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


def evidence_frame(rows: pd.DataFrame) -> pd.DataFrame:
    """Return evidence rows in the fixed column order, with missing columns empty.

    Parameters
    ----------
    rows : pd.DataFrame
        Evidence rows of one level.

    Returns
    -------
    pd.DataFrame
        ``rows`` with exactly ``EVIDENCE_COLUMNS``.

    Raises
    ------
    ValueError
        If ``rows`` has columns outside ``EVIDENCE_COLUMNS``.
    """
    unknown = [c for c in rows.columns if c not in EVIDENCE_COLUMNS]
    if unknown:
        raise ValueError(f"Unknown evidence columns: {unknown}")
    return rows.reindex(columns=EVIDENCE_COLUMNS)


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
    evidence_frame(rows).to_csv(path, sep="\t", index=False)
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


@contextlib.contextmanager
def staged_dir(out_dir: str) -> Iterator[str]:
    """Yield an empty directory that replaces ``out_dir`` only if the block succeeds.

    The staging directory is created beside ``out_dir``. On success the old
    ``out_dir`` is kept aside until the new one is in place; on any failure the
    staging directory is removed and ``out_dir`` is left as it was.

    Parameters
    ----------
    out_dir : str
        Final output directory; its parent must exist.

    Yields
    ------
    str
        Path of the staging directory.
    """
    parent, name = os.path.split(out_dir)
    staging = tempfile.mkdtemp(prefix=f".{name}_", dir=parent)
    # mkdtemp makes the directory private; give it the usual permissions
    umask = os.umask(0)
    os.umask(umask)
    os.chmod(staging, 0o777 & ~umask)
    try:
        yield staging
        backup = staging + ".old"
        if os.path.exists(out_dir):
            os.rename(out_dir, backup)
        try:
            os.replace(staging, out_dir)
        except OSError:
            if os.path.exists(backup):
                os.rename(backup, out_dir)
            raise
        shutil.rmtree(backup, ignore_errors=True)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
