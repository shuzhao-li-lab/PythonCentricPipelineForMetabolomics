"""
Blank filtering across separate Asari runs.

Feature IDs from different Asari runs are not compatible, so each feature in the
sample table is matched to the closest blank feature by m/z (ppm) and retention time.
A systematic retention time offset between the two runs is estimated first from
unambiguous m/z matches and applied before matching.

A feature is flagged when  blank_mean * ratio >= sample_mean,
where means are over non-zero values (same convention as pcpfm v1).
"""

import numpy as np
import pandas as pd

BLANK_COLUMNS = [
    "blank_match", "blank_mz", "blank_rtime", "blank_ppm", "blank_drt",
    "blank_mean", "sample_mean", "blank_ratio", "blank_flag",
]


def nonzero_mean(df):
    return df.where(df > 0).mean(axis=1).fillna(0.0)


def _candidates(mz_s, mz_b_sorted, ppm):
    """Index ranges into sorted blank m/z within ppm of each sample m/z."""
    tol = mz_s * ppm * 1e-6
    return np.searchsorted(mz_b_sorted, mz_s - tol, "left"), np.searchsorted(mz_b_sorted, mz_s + tol, "right")


def estimate_rt_offset(mz_s, rt_s, mz_b, rt_b, ppm=5, max_shift=60, min_pairs=20):
    """
    Median retention time difference (blank - sample) over features that have
    exactly one blank candidate within ppm and max_shift seconds. Returns (offset, n_pairs).
    """
    order = np.argsort(mz_b)
    mz_b, rt_b = mz_b[order], rt_b[order]
    lo, hi = _candidates(mz_s, mz_b, ppm)
    diffs = []
    for i in range(len(mz_s)):
        if hi[i] - lo[i] == 0:
            continue
        drt = rt_b[lo[i]:hi[i]] - rt_s[i]
        close = drt[np.abs(drt) <= max_shift]
        if len(close) == 1:
            diffs.append(close[0])
    if len(diffs) < min_pairs:
        return 0.0, len(diffs)
    return float(np.median(diffs)), len(diffs)


def match_features(sample_ft, blank_ft, ppm=5, rt_tol=10, rt_offset=0.0):
    """
    For each sample feature, the closest blank feature within ppm and rt_tol
    (after shifting blank rtime by -rt_offset). Distance is sqrt((dppm/ppm)^2 + (drt/rt_tol)^2).

    Returns a DataFrame indexed like sample_ft with blank_match, blank_mz, blank_rtime,
    blank_ppm, blank_drt (NaN / "" where there is no match).
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
        dppm = (mz_b[lo[i]:hi[i]] - mz_s[i]) / mz_s[i] * 1e6
        drt = rt_b[lo[i]:hi[i]] - rt_s[i]
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


def flag_blank_features(sample_ft, sample_cols, blank_values, ratio=3):
    """
    Add blank_mean, sample_mean, blank_ratio, blank_flag to sample_ft (a copy is returned).
    blank_values: Series aligned to sample_ft.index with the blank intensity for each feature (0 if none).
    """
    ft = sample_ft.copy()
    ft["blank_mean"] = blank_values.fillna(0.0).round(1)
    ft["sample_mean"] = nonzero_mean(ft[sample_cols]).round(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = ft["sample_mean"] / ft["blank_mean"]
    ft["blank_ratio"] = r.replace(np.inf, np.nan).round(2)
    ft["blank_flag"] = (ft["blank_mean"] > 0) & (ft["blank_mean"] * ratio >= ft["sample_mean"])
    return ft


def blank_filter_separate_runs(sample_ft, sample_cols, blank_ft, blank_cols,
                               ratio=3, ppm=5, rt_tol=10, estimate_offset=True):
    """Blank filtering between a sample-run table and a blank-run table. Returns (table, summary)."""
    offset, n_pairs = 0.0, 0
    if estimate_offset:
        offset, n_pairs = estimate_rt_offset(
            sample_ft["mz"].to_numpy(float), sample_ft["rtime"].to_numpy(float),
            blank_ft["mz"].to_numpy(float), blank_ft["rtime"].to_numpy(float), ppm=ppm,
        )
    matches = match_features(sample_ft, blank_ft, ppm=ppm, rt_tol=rt_tol, rt_offset=offset)
    blank_means = nonzero_mean(blank_ft[blank_cols])
    blank_means.index = blank_ft["id_number"]
    values = matches["blank_match"].map(blank_means).astype(float)

    ft = pd.concat([sample_ft, matches], axis=1)
    ft = flag_blank_features(ft, sample_cols, values, ratio)
    summary = {
        "mode": "separate_runs",
        "rt_offset_seconds": round(offset, 2),
        "rt_offset_pairs": n_pairs,
        "n_features": int(len(ft)),
        "n_matched": int((ft["blank_match"] != "").sum()),
        "n_flagged": int(ft["blank_flag"].sum()),
        "ratio": ratio, "ppm": ppm, "rt_tol": rt_tol,
    }
    return ft, summary


def blank_filter_same_run(ft, sample_cols, blank_cols, ratio=3):
    """Blank filtering when blanks and samples come from the same Asari run (same feature IDs)."""
    ft = ft.copy()
    ft["blank_match"] = np.where(nonzero_mean(ft[blank_cols]) > 0, ft["id_number"], "")
    for c in ("blank_mz", "blank_rtime", "blank_ppm", "blank_drt"):
        ft[c] = np.nan
    ft = flag_blank_features(ft, sample_cols, nonzero_mean(ft[blank_cols]), ratio)
    summary = {
        "mode": "same_run",
        "n_features": int(len(ft)),
        "n_matched": int((ft["blank_match"] != "").sum()),
        "n_flagged": int(ft["blank_flag"].sum()),
        "ratio": ratio,
    }
    return ft, summary
