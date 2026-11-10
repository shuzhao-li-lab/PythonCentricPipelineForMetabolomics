"""
Asari GC annotation for ``pcpfm annotate gc``.

Asari matches pseudospectra (co-eluting, correlated features) to an EI library
within a retention-index window computed from a Kovats file. Its table reader
takes columns by position, so the input file is written here in Asari's column
order from the named pcpfm columns.
"""

import contextlib
import os
import tempfile
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from ..tables import FEATURE_COLUMNS

# Asari parses these with int(), which rejects "12.0"
_INT_COLUMNS = ("parent_masstrack_id", "detection_counts")


def write_asari_input(
    table: pd.DataFrame, sample_cols: Sequence[str], path: str
) -> None:
    """Write a feature table in the column order Asari's GC reader expects.

    Parameters
    ----------
    table : pd.DataFrame
        pcpfm feature table with named Asari feature columns.
    sample_cols : Sequence[str]
        Intensity columns, written after the feature columns.
    path : str
        Output TSV path.

    Raises
    ------
    ValueError
        If a feature column is missing, has missing values, or an integer
        column holds non-integer values.

    Notes
    -----
    Columns that are neither feature nor sample columns (for example the blank
    filter columns) are left out, since Asari would read them as samples.
    """
    missing = [c for c in FEATURE_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(
            f"The table lacks Asari feature columns needed for GC annotation: {missing}"
        )
    clash = [c for c in sample_cols if c in FEATURE_COLUMNS]
    if clash:
        raise ValueError(f"Sample columns named like feature columns: {clash}")
    out = table[FEATURE_COLUMNS + list(sample_cols)].copy()
    # float turns pd.NA of nullable dtypes into NaN, which isfinite rejects
    numeric = (
        out[FEATURE_COLUMNS[1:]].apply(pd.to_numeric, errors="coerce").astype(float)
    )
    incomplete = [c for c in numeric.columns if not np.isfinite(numeric[c]).all()]
    if incomplete or out["id_number"].isna().any():
        raise ValueError(
            f"Feature columns with missing or non-finite values: {incomplete or ['id_number']}"
        )
    for c in _INT_COLUMNS:
        if not (out[c] == out[c].round()).all():
            raise ValueError(f"Column '{c}' must hold integers.")
        out[c] = out[c].astype(int)
    out.to_csv(path, sep="\t", index=False)


ASARI_LOG = "asari_gc_annotation.log"
# Asari names its outputs after this job name, e.g. Features_gc_annotation.tsv
JOB = "gc_annotation"
FEATURES_FILE = f"Features_{JOB}.tsv"
# fixed Asari settings: HR-EI fragment tolerance, RI spread of one compound, no PDFs
FIXED_PARAMS = {
    "ms2_tolerance_in_da": 0.005,
    "max_ri_delta": 100,
    "do_mirror_plot": False,
}


def run_asari_gc(
    table: pd.DataFrame,
    sample_cols: Sequence[str],
    kovats: str,
    db: str,
    params: dict[str, Any],
    out_dir: str,
) -> dict[str, Any]:
    """Annotate a GC feature table with Asari's EI library and Kovats search.

    Parameters
    ----------
    table : pd.DataFrame
        pcpfm feature table with named Asari feature columns.
    sample_cols : Sequence[str]
        Intensity columns of ``table``.
    kovats : str
        Kovats TSV file (RI, RT in minutes).
    db : str
        EI library in MSP or JSON format.
    params : dict[str, Any]
        Asari matching parameters: ``ri_tolerance``, ``score_cutoff_cosine``,
        ``score_cutoff_entropy``, ``corr_cutoff``, ``denovo`` (bool) and
        ``max_core_features``.
    out_dir : str
        Directory receiving Asari's outputs and its console log.

    Returns
    -------
    dict[str, Any]
        Path of Asari's Features table and log, the Asari version and every
        parameter passed to Asari.
    """
    import asari
    from asari.annotate import annotate_gcms_full

    passed = {**params, **FIXED_PARAMS}
    log_path = os.path.join(out_dir, ASARI_LOG)
    with tempfile.TemporaryDirectory() as tmp:
        infile = os.path.join(tmp, "asari_input.tsv")
        write_asari_input(table, sample_cols, infile)
        with (
            open(log_path, "w", encoding="utf-8") as log,
            contextlib.redirect_stdout(log),
            contextlib.redirect_stderr(log),
        ):
            annotate_gcms_full(
                infile, out_dir, kovats, db, project_name_handle=JOB, **passed
            )
    return {
        "features_tsv": os.path.join(out_dir, FEATURES_FILE),
        "log": log_path,
        "asari_version": asari.__version__,
        "params": passed,
    }
