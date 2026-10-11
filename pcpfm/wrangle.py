"""
Data wrangling from the filtered feature table to the analysis table:

    filtered table -> normalization -> batch correction -> imputation -> log2 -> analysis table

All functions take X as a features x samples DataFrame of intensities and return a new DataFrame.
``batches`` is a Series mapping sample -> batch label, or None.
"""

from typing import List, Literal, Tuple

import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_norm

STATS = {"median": np.median, "mean": np.mean, "min": np.min}


def _groups(X, batches):
    if batches is None:
        return {None: list(X.columns)}
    b = batches.reindex(X.columns).fillna("")
    return {k: list(v.index) for k, v in b.groupby(b)}


def normalize(X, method="none", batches=None, tic_percentile=0.9, mode="median"):
    """
    Sample normalization.

    tic: ported from v1 FeatureTable.TIC_normalize. Per-sample TIC is summed over features
        present in more than ``tic_percentile`` of samples; each sample is scaled to the
        median (or mean) TIC. With batches, normalize within each batch, then scale
        batch-level TICs to each other.
    median: scale each sample to the median of per-sample median non-zero intensities.

    TODO: more methods, e.g. probabilistic quotient normalization (PQN), internal standards,
    QC-based drift correction (LOESS on injection order).
    """
    if method == "none":
        return X
    stat = STATS[mode]
    if method == "median":
        med = X.where(X > 0).median(axis=0)
        return X / med * stat(med.to_numpy())
    if method != "tic":
        raise ValueError(f"Unknown normalization method: {method}")

    X = X.copy()
    batch_tics = {}
    for b, cols in _groups(X, batches).items():
        present = (X[cols] > 0).mean(axis=1) > tic_percentile
        tics = X.loc[present, cols].sum(axis=0)
        target = stat(tics.to_numpy())
        X[cols] = X[cols] * (target / tics)
        batch_tics[b] = target
    if len(batch_tics) > 1:
        target = stat(np.array(list(batch_tics.values())))
        for b, cols in _groups(X, batches).items():
            X[cols] = X[cols] * (target / batch_tics[b])
    return X


def batch_correct(
    X: pd.DataFrame,
    *,
    method: Literal["none", "combat"] = "none",
    batches: pd.Series | None = None,
) -> pd.DataFrame:
    """Correct batch effects in selected feature-table sample columns.

    Apply InMoose ComBat and floor corrected intensities at 1. Leave the
    input table unchanged, preserving metadata columns and row/column order.

    Parameters
    ----------
    X : pandas.DataFrame
        Feature table with one row per feature, and sample intensity columns.
    method : {"none", "combat"}, optional
        Correction method. The default, "none", returns X without correction.
        This is dictated by choices in CLI subparser. Currently supports Combat
    batches : pandas.Series or None, optional
        Batch labels indexed by sample column name, required for "combat".
        Labels are aligned to sample_columns regardless of their input order.
        Every selected sample must have a non-empty label, with at least two
        distinct batches represented. Extra sample labels are ignored.

    Returns
    -------
    pandas.DataFrame
        A corrected copy of X for "combat", or X itself for "none".
        Only the selected sample columns are modified in the copy.

    Raises
    ------
    ValueError
        If sample selection, table column names, the method, or batch labels
        are invalid; if the selected intensity table is empty, nonnumeric,
        or non-finite; or if ComBat produces non-finite values.

    Examples
    --------
    Correct four sample columns while retaining the feature metadata:

    >>> corrected = batch_correct(
    ...     feature_table,
    ...     method="combat",
    ...     batches=pd.Series({"S1": "A", "S2": "A", "S3": "B", "S4": "B"}),
    ... )
    """

    if method == "none":
        return X
    if method not in ["combat"]:  # Support more choices
        raise ValueError(f"Unknown batch correction method: {method}")
    if batches is None:
        raise ValueError(
            "Batch correction needs a batch column with at least two batches."
        )
    if not batches.index.is_unique:
        raise ValueError("Batch labels must have unique sample names.")
    labels = batches.reindex(X.columns)  # Make sure that the labels align
    if labels.isna().any() or labels.astype(str).str.strip().eq("").any():
        raise ValueError(
            "Batch correction needs a non-empty batch label for every sample."
        )
    if labels.nunique() < 2:
        raise ValueError(
            "Batch correction needs a batch column with at least two batches."
        )

    corrected = pycombat_norm(X.copy(), pd.factorize(labels)[0])
    if not np.isfinite(corrected.to_numpy()).all():
        raise ValueError(
            "ComBat produced non-finite intensities; check feature variation within batches."
        )
    return corrected


def impute(X, method="half-min", ratio=0.5, batches=None):
    """
    Fill zeros with ratio x the minimum non-zero value of each feature (v1 impute_missing_features),
    optionally per batch. Features never detected stay 0.
    """
    if method == "none":
        return X
    if method != "half-min":
        raise ValueError(f"Unknown imputation method: {method}")
    X = X.copy()
    for _, cols in _groups(X, batches).items():
        sub = X[cols]
        fill = sub.where(sub > 0).min(axis=1) * ratio
        X[cols] = sub.where(sub > 0, fill, axis=0).fillna(0)
    return X


def log2(X):
    return np.log2(X + 1)


def prepare_analysis_table(
    X: pd.DataFrame, meta, args
) -> Tuple[pd.DataFrame, List[str]]:
    """Apply wrangling steps per args. Returns (X, log lines).
    The input is expected to be a dataframe with only the sample columns.

    """
    batches = None
    if args.batch_field:
        if args.batch_field not in meta.columns:
            raise ValueError(
                f"Batch field '{args.batch_field}' is not in the metadata."
            )
        batches = meta[args.batch_field].reindex(X.columns)
    log = []
    X = normalize(X, args.normalize, batches)
    log.append(
        f"normalization: {args.normalize}"
        + (f" by {args.batch_field}" if batches is not None else "")
    )
    X = batch_correct(X, method=args.batch_correct, batches=batches)
    log.append(f"batch correction: {args.batch_correct}")
    X = impute(X, args.impute, batches=batches)
    log.append(f"imputation: {args.impute}")
    if not args.no_log2:
        X = log2(X)
        log.append("log2(x + 1) transform")
    return X, log
