"""
Blank filtering across separate Asari runs.

Feature IDs from different Asari runs are not compatible, so each feature in the
sample table is matched to the closest blank feature by m/z (ppm) and retention time.
A systematic retention time offset between the two runs is estimated first from
unambiguous m/z matches and applied before matching.

A feature is flagged when  blank_mean * ratio >= sample_mean,
where means are over non-zero values (same convention as pcpfm v1).
"""

from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import NDArray


def nonzero_mean(df: pd.DataFrame) -> pd.Series:
    """Compute each feature mean using only positive intensities.

    Parameters
    ----------
    df : pandas.DataFrame
        Feature-by-sample intensity table.

    Returns
    -------
    pd.Series
        Row means, with zero for rows without positive values.
    """
    return df.where(df > 0).mean(axis=1).fillna(0.0)


def _candidates(
    mz_s: NDArray[np.float64], mz_b_sorted: NDArray[np.float64], ppm: float
) -> tuple[NDArray[np.intp], NDArray[np.intp]]:
    """Find mass-window bounds in sorted blank feature masses.

    Parameters
    ----------
    mz_s : NDArray[np.float64]
        Sample feature masses.
    mz_b_sorted : NDArray[np.float64]
        Blank feature masses sorted in ascending order.
    ppm : float
        Symmetric mass tolerance in parts per million.

    Returns
    -------
    tuple[NDArray[np.intp], NDArray[np.intp]]
        Inclusive start and exclusive stop indices for each sample feature.
    """
    tol = mz_s * ppm * 1e-6
    return np.searchsorted(mz_b_sorted, mz_s - tol, "left"), np.searchsorted(
        mz_b_sorted, mz_s + tol, "right"
    )


def estimate_rt_offset(
    mz_s: NDArray[np.float64],
    rt_s: NDArray[np.float64],
    mz_b: NDArray[np.float64],
    rt_b: NDArray[np.float64],
    ppm: float = 5,
    max_shift: float = 60,
    min_pairs: int = 20,
) -> tuple[float, int]:
    """Estimate a global blank-minus-sample retention-time offset.

    Parameters
    ----------
    mz_s : NDArray[np.float64]
        Sample feature masses.
    rt_s : NDArray[np.float64]
        Sample retention times in seconds.
    mz_b : NDArray[np.float64]
        Blank feature masses.
    rt_b : NDArray[np.float64]
        Blank retention times in seconds.
    ppm : float, optional
        Mass tolerance in parts per million.
    max_shift : float, optional
        Maximum absolute retention-time difference in seconds.
    min_pairs : int, optional
        Minimum number of matched pairs required for an offset estimate.

    Returns
    -------
    tuple[float, int]
        Median offset in seconds and count of usable pairs; offset is zero below
        min_pairs.

    Notes
    -----
    Only unambiguous mass matches within max_shift contribute.
    """
    order = np.argsort(mz_b)
    mz_b, rt_b = mz_b[order], rt_b[order]
    lo, hi = _candidates(mz_s, mz_b, ppm)
    diffs = []
    for i in range(len(mz_s)):
        if hi[i] - lo[i] == 0:
            continue
        drt = rt_b[lo[i] : hi[i]] - rt_s[i]
        close = drt[np.abs(drt) <= max_shift]
        if len(close) == 1:
            diffs.append(close[0])
    if len(diffs) < min_pairs:
        return 0.0, len(diffs)
    return float(np.median(diffs)), len(diffs)


def match_features(
    sample_ft: pd.DataFrame,
    blank_ft: pd.DataFrame,
    ppm: float = 5,
    rt_tol: float = 10,
    rt_offset: float = 0.0,
) -> pd.DataFrame:
    """Match each sample feature to its closest eligible blank feature.

    Parameters
    ----------
    sample_ft : pandas.DataFrame
        Sample-run feature table with id_number, mz and rtime annotations.
    blank_ft : pandas.DataFrame
        Blank-run feature table with id_number, mz and rtime annotations.
    ppm : float, optional
        Mass tolerance in parts per million.
    rt_tol : float, optional
        Retention-time tolerance in seconds.
    rt_offset : float, optional
        Blank-minus-sample retention-time offset in seconds.

    Returns
    -------
    pd.DataFrame
        Blank match ID, mass, retention time, ppm error and time difference,
        indexed like sample_ft.

    Notes
    -----
    Candidates must satisfy both mass and adjusted retention-time tolerances.
    The smallest normalized squared distance wins. Unmatched rows have an empty
    ID and missing numeric fields.
    """
    mz_s = sample_ft["mz"].to_numpy(float)
    rt_s = sample_ft["rtime"].to_numpy(float)
    order = np.argsort(blank_ft["mz"].to_numpy(float))
    mz_b = blank_ft["mz"].to_numpy(float)[order]
    rt_b = blank_ft["rtime"].to_numpy(float)[order] - rt_offset
    ids_b = blank_ft["id_number"].to_numpy()[order]

    lo, hi = _candidates(mz_s, mz_b, ppm)
    best = np.full(len(mz_s), -1)
    for i in range(len(mz_s)):
        if hi[i] == lo[i]:
            continue
        dppm = (mz_b[lo[i] : hi[i]] - mz_s[i]) / mz_s[i] * 1e6
        drt = rt_b[lo[i] : hi[i]] - rt_s[i]
        ok = np.abs(drt) <= rt_tol
        if not ok.any():
            continue
        dist = np.where(ok, (dppm / ppm) ** 2 + (drt / rt_tol) ** 2, np.inf)
        best[i] = lo[i] + int(np.argmin(dist))

    has = best >= 0
    out = pd.DataFrame(index=sample_ft.index)
    out["blank_match"] = np.where(has, ids_b[best], "")
    out["blank_mz"] = np.where(has, mz_b[best], np.nan)
    out["blank_rtime"] = np.where(has, rt_b[best] + rt_offset, np.nan)
    out["blank_ppm"] = np.where(has, (mz_b[best] - mz_s) / mz_s * 1e6, np.nan).round(2)
    out["blank_drt"] = np.where(has, rt_b[best] - rt_s, np.nan).round(2)
    return out


def flag_blank_features(
    sample_ft: pd.DataFrame,
    sample_cols: list[str],
    blank_values: pd.Series,
    ratio: float = 3,
) -> pd.DataFrame:
    """Flag features whose blank signal is large relative to sample signal.

    Parameters
    ----------
    sample_ft : pandas.DataFrame
        Sample-run feature table with id_number, mz and rtime annotations.
    sample_cols : list[str]
        Sample intensity column names in the feature table, in the desired
        order.
    blank_values : pandas.Series
        Blank feature means aligned to the sample feature index.
    ratio : float, optional
        Sample-to-blank threshold for flagging background features.

    Returns
    -------
    pd.DataFrame
        Copy with sample_mean, blank_mean, blank_ratio and blank_flag
        annotations.

    Notes
    -----
    A feature is flagged when its blank mean is positive and blank_mean * ratio
    >= sample_mean. Comparisons use unrounded means.
    """
    ft = sample_ft.copy()
    blank_mean = blank_values.fillna(0.0)
    sample_mean = nonzero_mean(ft[sample_cols])
    with np.errstate(divide="ignore", invalid="ignore"):
        r = sample_mean / blank_mean
    # rounded for readable tables only; the ratio and flag use the unrounded means
    ft["blank_mean"] = blank_mean.round(1)
    ft["sample_mean"] = sample_mean.round(1)
    ft["blank_ratio"] = r.replace(np.inf, np.nan).round(2)
    ft["blank_flag"] = (blank_mean > 0) & (blank_mean * ratio >= sample_mean)
    return ft


def blank_filter_separate_runs(
    sample_ft: pd.DataFrame,
    sample_cols: list[str],
    blank_ft: pd.DataFrame,
    blank_cols: list[str],
    ratio: float = 3,
    ppm: float = 5,
    rt_tol: float = 10,
    estimate_offset: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Match separate sample and blank runs and annotate background features.

    Parameters
    ----------
    sample_ft : pandas.DataFrame
        Sample-run feature table with id_number, mz and rtime annotations.
    sample_cols : list[str]
        Sample intensity column names in the feature table, in the desired
        order.
    blank_ft : pandas.DataFrame
        Blank-run feature table with id_number, mz and rtime annotations.
    blank_cols : list[str]
        Blank intensity column names in the feature table.
    ratio : float, optional
        Sample-to-blank threshold for flagging background features.
    ppm : float, optional
        Mass tolerance in parts per million.
    rt_tol : float, optional
        Retention-time tolerance in seconds.
    estimate_offset : bool, optional
        Whether to estimate a global retention-time offset before matching.

    Returns
    -------
    tuple[pd.DataFrame, dict[str, Any]]
        Annotated sample table and matching/filtering summary. No feature rows
        are removed.
    """
    offset, n_pairs = 0.0, 0
    if estimate_offset:
        offset, n_pairs = estimate_rt_offset(
            sample_ft["mz"].to_numpy(float),
            sample_ft["rtime"].to_numpy(float),
            blank_ft["mz"].to_numpy(float),
            blank_ft["rtime"].to_numpy(float),
            ppm=ppm,
        )
    matches = match_features(
        sample_ft, blank_ft, ppm=ppm, rt_tol=rt_tol, rt_offset=offset
    )
    blank_means = nonzero_mean(blank_ft[blank_cols])
    blank_means.index = blank_ft["id_number"]
    values = matches["blank_match"].map(blank_means).astype(float)

    ft = pd.concat([sample_ft, matches], axis=1)
    ft = flag_blank_features(ft, sample_cols, values, ratio)
    summary = {
        "mode": "separate_runs",
        "rt_offset_seconds": round(offset, 2),
        "rt_offset_pairs": n_pairs,
        "n_features": len(ft),
        "n_matched": int((ft["blank_match"] != "").sum()),
        "n_flagged": int(ft["blank_flag"].sum()),
        "ratio": ratio,
        "ppm": ppm,
        "rt_tol": rt_tol,
    }
    return ft, summary


def blank_filter_same_run(
    ft: pd.DataFrame, sample_cols: list[str], blank_cols: list[str], ratio: float = 3
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Annotate background features using blank columns in the same table.

    Parameters
    ----------
    ft : pandas.DataFrame
        Feature table containing sample intensities and available feature
        annotations.
    sample_cols : list[str]
        Sample intensity column names in the feature table, in the desired
        order.
    blank_cols : list[str]
        Blank intensity column names in the feature table.
    ratio : float, optional
        Sample-to-blank threshold for flagging background features.

    Returns
    -------
    tuple[pd.DataFrame, dict[str, Any]]
        Annotated feature table and blank-filter summary. No feature rows are
        removed.
    """
    ft = ft.copy()
    ft["blank_match"] = np.where(nonzero_mean(ft[blank_cols]) > 0, ft["id_number"], "")
    for c in ("blank_mz", "blank_rtime", "blank_ppm", "blank_drt"):
        ft[c] = np.nan
    ft = flag_blank_features(ft, sample_cols, nonzero_mean(ft[blank_cols]), ratio)
    summary = {
        "mode": "same_run",
        "n_features": len(ft),
        "n_matched": int((ft["blank_match"] != "").sum()),
        "n_flagged": int(ft["blank_flag"].sum()),
        "ratio": ratio,
    }
    return ft, summary
