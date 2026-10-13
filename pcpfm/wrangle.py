"""
Transform filtered feature intensities into an analysis table.

    filtered table -> normalization -> batch correction -> imputation -> log2 -> analysis table

All transformations take X as a features x samples DataFrame of intensities.
Inputs are not modified; disabled transformations return X itself.
``batches`` is a Series mapping sample -> batch label, or None.
"""

from argparse import Namespace
from collections.abc import Hashable
from typing import Literal

import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_norm

STATS = {"median": np.median, "mean": np.mean, "min": np.min}


def _groups(
    X: pd.DataFrame, batches: pd.Series | None
) -> dict[Hashable | None, list[str]]:
    """Group intensity columns by their aligned batch labels.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    batches : pandas.Series | None
        Batch labels indexed by sample column name, or None for a single group.

    Returns
    -------
    dict[Hashable | None, list[str]]
        Batch labels mapped to sample columns; one None group when batches are
        absent.

    Notes
    -----
    Missing batch labels form the empty-string group.
    """
    if batches is None:
        return {None: list(X.columns)}
    b = batches.reindex(X.columns).fillna("")
    return {k: list(v.index) for k, v in b.groupby(b)}


def normalize(
    X: pd.DataFrame,
    method: Literal["none", "median", "tic"] = "none",
    batches: pd.Series | None = None,
    tic_percentile: float = 0.9,
    mode: Literal["median", "mean", "min"] = "median",
) -> pd.DataFrame:
    """Normalize sample intensities by median intensity or shared-feature TIC.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    method : Literal['none', 'median', 'tic'], optional
        Normalization method; "none" skips normalization.
    batches : pandas.Series | None, optional
        Batch labels indexed by sample column name, or None for a single group.
    tic_percentile : float, optional
        Minimum fraction of samples detecting a feature for inclusion in TIC.
    mode : Literal['median', 'mean', 'min'], optional
        Statistic used to choose the target intensity or TIC.

    Returns
    -------
    pd.DataFrame
        Normalized intensities, or X itself for method="none". The input is not
        modified.

    Raises
    ------
    ValueError
        If the method is unknown or a sample cannot be scaled because it has no detected features, no shared TIC features, or zero TIC.

    Notes
    -----
    TIC normalization sums features detected in strictly more than
    tic_percentile of a batch, scales within batches, then scales batch TICs to
    a common target. Median normalization uses positive per-sample medians
    across all samples. Zeros remain zero.
    """
    if method not in ["none", "median", "tic"]:
        raise ValueError(f"Unknown normalization method: {method}")

    if method == "none":
        return X
    stat = STATS[mode]
    if method == "median":
        med = X.where(X > 0).median(axis=0)
        # a sample with no detected features has a NaN median, which would make the
        # target (median of medians) NaN for every sample
        empty = med.index[med.isna()].tolist()
        if empty:
            raise ValueError(
                f"median normalization: samples with no detected features: {empty}"
            )
        return X / med * stat(med.to_numpy())

    X = X.copy()
    batch_tics = {}
    for b, cols in _groups(X, batches).items():
        present = (X[cols] > 0).mean(axis=1) > tic_percentile
        # without common features there is no TIC to scale to
        if not present.any():
            where = "all samples" if b is None else f"samples in batch '{b}'"
            raise ValueError(
                f"TIC normalization: no feature is present in more than {tic_percentile:.0%} "
                f"of {where}; use --normalize median or none."
            )
        tics = X.loc[present, cols].sum(axis=0)
        # a sample with TIC 0 would be divided by zero
        zero = tics.index[tics == 0].tolist()
        if zero:
            raise ValueError(
                f"TIC normalization: samples with a TIC of 0 cannot be scaled: {zero}"
            )
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
    """Correct batch effects in an intensity-only feature matrix using ComBat.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    method : Literal['none', 'combat'], optional
        Correction method, also constrained by the CLI choices.
    batches : pandas.Series | None, optional
        Batch labels indexed by sample column name. Required for "combat"; every
        sample needs a non-empty label.

    Returns
    -------
    pd.DataFrame
        ComBat-corrected intensities with the input labels, or X itself for
        method="none".

    Raises
    ------
    ValueError
        If the method is unknown, the ComBat input is empty, batch labels are missing or duplicated, fewer than two batches remain, or ComBat returns non-finite values.

    Notes
    -----
    Pass only sample intensity columns, excluding feature annotations. Batch
    labels are aligned to X.columns; extra labels are ignored. The input is not
    modified. Corrected values are returned without clipping. InMoose validates
    the numeric input and may raise additional errors.
    """
    if method == "none":
        return X
    if method not in ["combat"]:  # Support more choices
        raise ValueError(f"Unknown batch correction method: {method}")
    if X.empty:
        raise ValueError("Batch correction input table is empty")
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


def impute(
    X: pd.DataFrame,
    method: Literal["none", "half-min"] = "half-min",
    ratio: float = 0.5,
    batches: pd.Series | None = None,
) -> pd.DataFrame:
    """Replace undetected intensities with a fraction of each feature minimum.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    method : Literal['none', 'half-min'], optional
        Imputation method; "none" leaves the data unchanged.
    ratio : float, optional
        Fraction of the minimum positive intensity used for filling.
    batches : pandas.Series | None, optional
        Batch labels indexed by sample column name, or None for a single group.

    Returns
    -------
    pd.DataFrame
        Imputed copy, or X itself for method="none".

    Raises
    ------
    ValueError
        If the imputation method is unknown.

    Notes
    -----
    For each batch, nonpositive and missing values are replaced by ratio times
    the minimum positive feature intensity. Features never detected in a batch
    remain zero.
    """
    if method == "none":
        return X
    if method != "half-min":
        raise ValueError(f"Unknown imputation method: {method}")
    X = X.copy()
    for cols in _groups(X, batches).values():
        sub = X[cols]
        fill = sub.where(sub > 0).min(axis=1) * ratio
        X[cols] = sub.where(sub > 0, fill, axis=0).fillna(0)
    return X


def log2(X: pd.DataFrame) -> pd.DataFrame:
    """Vectorized Transform intensities with log2(x + 1).

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.

    Returns
    -------
    pd.DataFrame
        Transformed intensities with the same row and column labels.
    """
    return np.log2(X + 1)


def prepare_analysis_table(
    X: pd.DataFrame, meta: pd.DataFrame, args: Namespace
) -> tuple[pd.DataFrame, list[str]]:
    """Normalize, batch-correct, impute and optionally log-transform intensities.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    meta : pandas.DataFrame
        Sample metadata indexed by sample ID.
    args : Namespace
        Parsed command options for this pipeline step.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        Analysis intensity table and descriptions of the applied steps.

    Raises
    ------
    ValueError
        If the batch field is missing, a wrangling step rejects its input, or the final table contains non-finite values.
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
    if not np.isfinite(X.to_numpy(dtype=float)).all():
        raise ValueError(
            f"analysis table contains NaN or infinite values after: {'; '.join(log)}"
        )
    return X, log
