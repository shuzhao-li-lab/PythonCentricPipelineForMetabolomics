"""
Data wrangling from the filtered feature table to the analysis table:

    filtered table -> normalization -> batch correction -> imputation -> log2 -> analysis table

All functions take X as a features x samples DataFrame of intensities and return a new DataFrame.
``batches`` is a Series mapping sample -> batch label, or None.
"""

import numpy as np
import pandas as pd

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


def batch_correct(X, method="none", batches=None):
    """
    Batch correction. Placeholder.

    v1 used ComBat (pycombat) on the intensity table with batches from a metadata field,
    and set negative values afterwards to 1. To be implemented in v2, candidates:
    ComBat (e.g. inmoose.pycombat), QC-sample based correction (QC-RLSC / LOESS per batch),
    or WaveICA.
    """
    if method == "none":
        return X
    if batches is None or batches.reindex(X.columns).nunique() < 2:
        raise ValueError("Batch correction needs a batch column with at least two batches.")
    raise NotImplementedError(f"Batch correction '{method}' is not implemented yet in pcpfm 2.")


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


def prepare_analysis_table(X, meta, args):
    """Apply wrangling steps per args. Returns (X, log lines)."""
    batches = None
    if args.batch_field:
        if args.batch_field not in meta.columns:
            raise ValueError(f"Batch field '{args.batch_field}' is not in the metadata.")
        batches = meta[args.batch_field].reindex(X.columns)
    log = []
    X = normalize(X, args.normalize, batches)
    log.append(f"normalization: {args.normalize}" + (f" by {args.batch_field}" if batches is not None else ""))
    X = batch_correct(X, args.batch_correct, batches)
    log.append(f"batch correction: {args.batch_correct}")
    X = impute(X, args.impute, batches=batches)
    log.append(f"imputation: {args.impute}")
    if not args.no_log2:
        X = log2(X)
        log.append("log2(x + 1) transform")
    return X, log
