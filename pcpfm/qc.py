"""
``pcpfm qc``: quality control on a feature table, using metadata when available.

Per sample (sample_qc.tsv):
    detected features, missing rate, total and median log2 intensity,
    median Spearman correlation to other samples, PCA Hotelling T2,
    Asari RT-alignment landmark peaks; robust z-scores (median/MAD) and flags.
Per feature (feature_qc.tsv):
    Asari peak quality (snr, goodness_fitting, cSelectivity), detection rate,
    blank flag, CV and D-ratio from QC samples (e.g. pooled), dilution-series linearity.

Tables written:
    tables/feature_table_filtered.tsv   blank-flagged and undetected features removed; raw intensities
    tables/feature_table_analysis.tsv   filtered, then normalized, batch-corrected, imputed and log2
                                        transformed (see wrangle.py); entry point for `pcpfm analyze`
The QC report is qc/qc_report.html.
"""

from __future__ import annotations

import os
from argparse import Namespace
from collections.abc import Iterable
from typing import Any

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray
from scipy import stats

from .metadata import is_blank, is_study_sample
from .project import COMPLETED, FAILED, ProjectRecord
from .tables import load_table, parse_selector, sample_info, select

ZSCORE_LIMIT = 3.0
# metadata columns that define which samples are comparable replicates for QC flags
QC_STRATA = ("sample_type", "group")


def robust_z(x: pd.Series | ArrayLike) -> pd.Series:
    """Standardize values using the median and scaled median absolute deviation.

    Parameters
    ----------
    x : pandas.Series | ArrayLike
        One-dimensional numeric observations.

    Returns
    -------
    pd.Series
        Robust z-scores with the input index when x is a Series.

    Notes
    -----
    The MAD is scaled by 1.4826. Zero or undefined MAD yields zero for finite
    inputs; missing values remain missing.
    """
    x = pd.Series(x, dtype=float)
    mad = 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))
    if not mad or np.isnan(mad):
        return x * 0.0
    return (x - np.nanmedian(x)) / mad


# def stratified_robust_z(
#     values: pd.Series | ArrayLike, strata: Iterable[str], min_size: int = 3
# ) -> pd.Series:
#     """Compute robust z-scores within comparable sample strata.
#
#     Parameters
#     ----------
#     values : pandas.Series | ArrayLike
#         One-dimensional numeric observations.
#     strata : Iterable[str]
#         Stratum labels in the same order as values.
#     min_size : int, optional
#         Minimum number of observations in a stratum.
#
#     Returns
#     -------
#     pd.Series
#         Z-scores in input order; NaN for strata smaller than min_size.
#
#     Notes
#     -----
#     Strata are aligned positionally. Each stratum uses its own median and scaled
#     MAD; zero MAD gives zero for finite values.
#     """
#     values = pd.Series(values, dtype=float)
#     strata = pd.Series(list(strata), index=values.index)
#     z = values.groupby(strata).transform(robust_z)
#     return z.where(strata.map(strata.value_counts()) >= min_size)


def qc_strata(meta: pd.DataFrame) -> pd.Series:
    """Build sample strata from available sample_type and group columns.

    Parameters
    ----------
    meta : pandas.DataFrame
        Sample metadata indexed by sample ID.

    Returns
    -------
    pd.Series
        Labels indexed by sample ID, joining metadata values with "|".
    """
    cols = [c for c in QC_STRATA if c in meta.columns]
    return meta[cols].fillna("").astype(str).agg("|".join, axis=1)


def rowwise_pearson(
    Y: NDArray[np.float64], x: NDArray[np.float64], min_n: int = 4
) -> NDArray[np.float64]:
    """Correlate each feature row with a sample-level numeric vector.

    Parameters
    ----------
    Y : NDArray[np.float64]
        Feature-by-sample numeric array, with NaN for missing observations.
    x : NDArray[np.float64]
        Finite numeric predictor in sample-column order.
    min_n : int, optional
        Minimum number of observed intensities per feature.

    Returns
    -------
    NDArray[np.float64]
        Pearson coefficients; NaN for constant rows or rows with fewer than
        min_n observations.

    Notes
    -----
    NaN entries in Y are excluded separately for each row. The predictor x must
    contain finite values.
    """
    mask = ~np.isnan(Y)
    n = mask.sum(axis=1)
    Xb = np.where(mask, x, 0.0)
    Yz = np.where(mask, Y, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        xm = Xb.sum(axis=1) / n
        ym = Yz.sum(axis=1) / n
        dx = np.where(mask, Xb - xm[:, None], 0.0)
        dy = np.where(mask, Yz - ym[:, None], 0.0)
        r = (dx * dy).sum(axis=1) / np.sqrt((dx**2).sum(axis=1) * (dy**2).sum(axis=1))
    r[n < min_n] = np.nan
    return r


# ---------------- computations ----------------


def pca(
    X: pd.DataFrame, min_detection: float = 0.5, max_components: int = 10
) -> tuple[pd.DataFrame, NDArray[np.float64], pd.Series, float, int] | None:
    """Compute PCA and Hotelling statistics from detected feature intensities.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    min_detection : float, optional
        Minimum fraction of samples detecting a feature.
    max_components : int, optional
        Maximum number of principal components to retain.

    Returns
    -------
    tuple[pd.DataFrame, NDArray[np.float64], pd.Series, float, int] | None
        Sample scores, explained variance ratios, sample T-squared values, 95%
        T-squared limit and number of components used for T-squared. None when
        fewer than three retained features or four samples are available.

    Notes
    -----
    Retain features meeting min_detection, replace undetected values by half the
    positive feature minimum, then log2-transform and autoscale before SVD.
    """
    X = X[(X > 0).mean(axis=1) >= min_detection]
    n = X.shape[1]
    if len(X) < 3 or n < 4:
        return None
    half_min = X.where(X > 0).min(axis=1) / 2
    L = np.log2(X.where(X > 0, half_min, axis=0)).T.to_numpy()
    sd = L.std(axis=0, ddof=1)
    L = (L[:, sd > 0] - L[:, sd > 0].mean(axis=0)) / sd[sd > 0]
    U, S, _ = np.linalg.svd(L, full_matrices=False)
    ncomp = min(max_components, n - 1)
    var = S**2 / np.sum(S**2)
    scores = U[:, :ncomp] * S[:ncomp]
    # Hotelling T2 on the components explaining 80% of variance (2 to 5 components)
    k = int(np.clip(np.searchsorted(np.cumsum(var), 0.8) + 1, 2, min(5, ncomp)))
    lam = S[:k] ** 2 / (n - 1)
    t2 = np.sum(scores[:, :k] ** 2 / lam, axis=1)
    limit = k * (n - 1) / (n - k) * stats.f.ppf(0.95, k, n - k) if n > k else np.inf
    cols = [f"PC{i + 1}" for i in range(ncomp)]
    return (
        pd.DataFrame(scores, index=X.columns, columns=cols),
        var[:ncomp],
        pd.Series(t2, index=X.columns),
        float(limit),
        k,
    )


def sample_metrics(
    X: pd.DataFrame,
    meta: pd.DataFrame,
    alignment: dict[str, dict[str, Any]] | None = None,
    pca_result: tuple[pd.DataFrame, NDArray[np.float64], pd.Series, float, int]
    | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute sample intensity, correlation and outlier metrics.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    meta : pandas.DataFrame
        Sample metadata indexed by sample ID.
    alignment : dict[str, dict[str, Any]] | None, optional
        Optional Asari alignment metrics keyed by sample ID.
    pca_result : tuple[pandas.DataFrame, NDArray[np.float64], pandas.Series, float, int] | None, optional
        PCA scores, variance ratios, T-squared values, threshold and component
        count, or None when unavailable.

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Sample metrics indexed by sample ID, and the sample Spearman correlation
        matrix.
    """
    df = pd.DataFrame(index=X.columns)
    df.index.name = "sample"
    df["sample_type"] = meta["sample_type"]
    for c in ("group", "batch", "injection_order", "dilution_factor"):
        if c in meta.columns:
            df[c] = meta[c]
    detected = X > 0
    df["n_features"] = detected.sum(axis=0)
    df["missing_rate"] = (1 - detected.mean(axis=0)).round(4)
    df["total_intensity"] = X.sum(axis=0)
    df["median_log2_intensity"] = np.log2(X.where(detected)).median(axis=0).round(3)
    L = np.log2(X.where(detected))
    corr = L.corr(method="spearman", min_periods=20)
    corr = corr.mask(np.eye(len(corr), dtype=bool))
    df["median_correlation"] = corr.median(axis=1).round(4)

    if pca_result is not None:
        scores, _, t2, _limit, _ = pca_result
        df["PC1"], df["PC2"] = scores["PC1"].round(3), scores["PC2"].round(3)
        df["pca_T2"] = t2.round(2)
    if alignment:
        df["asari_landmark_peaks"] = [
            alignment.get(s, {}).get("landmark_peaks", np.nan) for s in df.index
        ]
        df["asari_rt_aligned"] = [
            alignment.get(s, {}).get("aligned", True) for s in df.index
        ]

    # compare each sample only with others of the same sample_type and group
    # strata = qc_strata(meta.reindex(df.index))
    # df["qc_stratum"] = strata
    df["z_n_features"] = robust_z(df["n_features"]).round(2)
    df["z_total_intensity"] = robust_z(np.log2(df["total_intensity"] + 1)).round(2)
    df["z_median_correlation"] = robust_z(df["median_correlation"]).round(2)

    # a NaN z-score (stratum too small) compares False, so it never raises a flag
    flags = {s: [] for s in df.index}
    for s in df.index:
        if df.at[s, "z_n_features"] < -ZSCORE_LIMIT:
            flags[s].append("low_feature_count")
        if abs(df.at[s, "z_total_intensity"]) > ZSCORE_LIMIT:
            flags[s].append("total_intensity_outlier")
        if df.at[s, "z_median_correlation"] < -ZSCORE_LIMIT:
            flags[s].append("low_correlation")
        if pca_result is not None and df.at[s, "pca_T2"] > pca_result[3]:
            flags[s].append("pca_outlier")
        if alignment and not df.at[s, "asari_rt_aligned"]:
            flags[s].append("rt_alignment_failed")
    df["flags"] = [";".join(flags[s]) for s in df.index]
    return df, corr


def feature_metrics(
    ft: pd.DataFrame,
    X: pd.DataFrame,
    meta: pd.DataFrame,
    qc_samples: list[str],
    study_samples: list[str],
) -> pd.DataFrame:
    """Compute detection, peak-quality and available replicate QC metrics.

    Parameters
    ----------
    ft : pandas.DataFrame
        Feature table containing sample intensities and available feature
        annotations.
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    meta : pandas.DataFrame
        Sample metadata indexed by sample ID.
    qc_samples : list[str]
        Sample IDs selected as technical QC replicates.
    study_samples : list[str]
        Study sample IDs used to estimate biological variation for D-ratio.

    Returns
    -------
    pd.DataFrame
        Feature metrics in feature-table order, including CV, D-ratio and
        dilution correlation when supported.
    """
    cols = [
        c
        for c in (
            "id_number",
            "mz",
            "rtime",
            "snr",
            "goodness_fitting",
            "cSelectivity",
            "blank_flag",
        )
        if c in ft.columns
    ]
    df = ft[cols].copy()
    detected = X > 0
    df["detection_rate"] = detected.mean(axis=1).round(4).to_numpy()
    df["mean_intensity"] = X.where(detected).mean(axis=1).round(1).to_numpy()

    if len(qc_samples) >= 3:
        Q = X[qc_samples].where(X[qc_samples] > 0)
        df["qc_detection_rate"] = (X[qc_samples] > 0).mean(axis=1).round(4).to_numpy()
        cv = (Q.std(axis=1, ddof=1) / Q.mean(axis=1) * 100).where(Q.count(axis=1) >= 2)
        df["qc_cv"] = cv.round(2).to_numpy()
        if len(study_samples) >= 3:
            St = X[study_samples].where(X[study_samples] > 0)
            d = Q.std(axis=1, ddof=1) / St.std(axis=1, ddof=1) * 100
            df["d_ratio"] = d.replace(np.inf, np.nan).round(2).to_numpy()

    if "dilution_factor" in meta.columns:
        dil = pd.to_numeric(meta["dilution_factor"], errors="coerce")
        dil = dil[dil > 0]
        if dil.nunique() >= 3:
            L = np.log2(X[dil.index].where(X[dil.index] > 0)).to_numpy()
            r = rowwise_pearson(L, np.log2(1 / dil.to_numpy(float)), min_n=4)
            df["dilution_r"] = np.round(r, 3)
    return df


def summarize(
    ft: pd.DataFrame,
    smet: pd.DataFrame,
    fmet: pd.DataFrame,
    qc_samples: list[str],
    pca_result: tuple[pd.DataFrame, NDArray[np.float64], pd.Series, float, int] | None,
    blank_summary: dict[str, Any] | None,
) -> dict[str, Any]:
    """Summarize sample and feature QC for the project record and report.

    Parameters
    ----------
    ft : pandas.DataFrame
        Feature table containing sample intensities and available feature
        annotations.
    smet : pandas.DataFrame
        Sample QC metrics indexed by sample ID.
    fmet : pandas.DataFrame
        Feature QC metrics in feature-table row order.
    qc_samples : list[str]
        Sample IDs selected as technical QC replicates.
    pca_result : tuple[pandas.DataFrame, NDArray[np.float64], pandas.Series, float, int] | None
        PCA scores, variance ratios, T-squared values, threshold and component
        count, or None when unavailable.
    blank_summary : dict[str, Any] | None
        Recorded blank-filter summary, if available.

    Returns
    -------
    dict[str, Any]
        Counts, flagged samples, available QC statistics and blank-filter
        information.
    """
    s = {
        "n_features": len(ft),
        "n_samples": len(smet),
        "sample_types": smet["sample_type"].value_counts().to_dict(),
        "median_features_per_sample": float(smet["n_features"].median()),
        "median_detection_rate": float(fmet["detection_rate"].median()),
        "flagged_samples": {k: v for k, v in smet["flags"].items() if v},
    }
    if "blank_flag" in fmet.columns:
        s["blank_flagged_features"] = int(fmet["blank_flag"].sum())
    if blank_summary:
        s["blank_filter"] = blank_summary
    if "qc_cv" in fmet.columns:
        cv = fmet["qc_cv"].dropna()
        s["qc_samples"] = list(qc_samples)
        s["qc_cv_median"] = round(float(cv.median()), 2) if len(cv) else None
        s["qc_cv_below_20_pct"] = (
            round(float((cv < 20).mean() * 100), 1) if len(cv) else None
        )
        s["qc_cv_below_30_pct"] = (
            round(float((cv < 30).mean() * 100), 1) if len(cv) else None
        )
    if "d_ratio" in fmet.columns:
        d = fmet["d_ratio"].dropna()
        s["d_ratio_below_50_pct"] = (
            round(float((d < 50).mean() * 100), 1) if len(d) else None
        )
    if "dilution_r" in fmet.columns:
        r = fmet["dilution_r"].dropna()
        s["dilution_features_tested"] = len(r)
        s["dilution_r_above_0.8_pct"] = (
            round(float((r > 0.8).mean() * 100), 1) if len(r) else None
        )
    if pca_result is not None:
        s["pca_explained_variance"] = [round(float(v), 4) for v in pca_result[1][:5]]
        s["pca_T2_limit_95"] = round(pca_result[3], 2)
        s["pca_T2_components"] = pca_result[4]
    return s


# ---------------- wrangling ----------------


def filter_table(
    ft: pd.DataFrame,
    sample_cols: list[str],
    fmet: pd.DataFrame,
    smet: pd.DataFrame,
    args: Namespace,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Filter features and optionally flagged samples without changing intensities.

    Parameters
    ----------
    ft : pandas.DataFrame
        Feature table containing sample intensities and available feature
        annotations.
    sample_cols : list[str]
        Sample intensity column names in the feature table, in the desired
        order.
    fmet : pandas.DataFrame
        Feature QC metrics in feature-table row order.
    smet : pandas.DataFrame
        Sample QC metrics indexed by sample ID.
    args : Namespace
        Parsed command options for this pipeline step.

    Returns
    -------
    tuple[pd.DataFrame, list[str], list[str]]
        Filtered feature table, remaining sample columns and filtering log.

    Notes
    -----
    Applies blank flags, detection rate and optional QC CV thresholds. Feature
    metrics must follow feature-table row order.
    """
    log = [f"start: {len(ft)} features, {len(sample_cols)} samples"]
    keep = pd.Series(True, index=ft.index)
    if not args.keep_blank_features and "blank_flag" in ft.columns:
        keep &= ~ft["blank_flag"].astype(bool)
        log.append(f"remove blank-flagged features: {int((~keep).sum())} removed")
    before = int(keep.sum())
    keep &= fmet["detection_rate"].to_numpy() >= args.min_detection
    keep &= fmet["detection_rate"].to_numpy() > 0
    log.append(
        f"detection rate >= {args.min_detection}: {before - int(keep.sum())} removed"
    )
    if args.max_qc_cv and "qc_cv" in fmet.columns:
        before = int(keep.sum())
        keep &= (fmet["qc_cv"].fillna(np.inf) <= args.max_qc_cv).to_numpy()
        log.append(f"QC CV <= {args.max_qc_cv}%: {before - int(keep.sum())} removed")
    out = ft[keep.to_numpy()].copy()

    if args.drop_flagged_samples:
        drop = [s for s in sample_cols if smet.at[s, "flags"]]
        sample_cols = [s for s in sample_cols if s not in drop]
        out = out.drop(columns=drop)
        log.append(f"drop flagged samples: {drop}")

    log.append(f"end: {len(out)} features, {len(sample_cols)} samples")
    return out, sample_cols, log


# ---------------- entry point ----------------


def run(args: Namespace, record: ProjectRecord) -> None:
    """Run QC and persist its completed or failed status.

    Parameters
    ----------
    args : Namespace
        Parsed command options for this pipeline step.
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.

    Notes
    -----
    The step is marked running before work starts. Exceptions are recorded as
    failures and re-raised.
    """
    params = {k: v for k, v in vars(args).items() if k not in ("func", "argv")}
    record.start_step("qc", params)
    try:
        summary = _qc(args, record)
    except Exception as e:
        record.finish_step("qc", FAILED, message=str(e))
        raise
    record.finish_step("qc", COMPLETED, summary=summary)


def _qc(args: Namespace, record: ProjectRecord) -> dict[str, Any]:
    """Execute QC and register generated tables and reports.

    Parameters
    ----------
    args : Namespace
        Parsed command options for this pipeline step.
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.

    Returns
    -------
    dict[str, Any]
        Step summary for the project record.
    """
    from . import report, wrangle

    ft, sample_cols = load_table(record, args.table)
    all_cols = list(sample_cols)
    X = ft[sample_cols].astype(float).fillna(0)
    meta = sample_info(record, sample_cols)
    blank_like = [s for s in sample_cols if is_blank(meta.at[s, "sample_type"])]
    if blank_like:
        print(f"[pcpfm] Excluding {len(blank_like)} blank samples from QC metrics.")
        sample_cols = [s for s in sample_cols if s not in blank_like]
        X, meta = X[sample_cols], meta.loc[sample_cols]

    qc_samples = select(meta, parse_selector(args.qc_samples))
    # study samples for the D-ratio denominator: blanks are already removed above;
    # QC samples are the numerator; other non-study types (standards, NIST...) and
    # dilution-series injections are designed variation, not biological variation
    dilution = (
        meta.get("dilution_factor", pd.Series("", index=meta.index))
        .fillna("")
        .astype(str)
    )
    study_samples = [
        s
        for s in sample_cols
        if s not in qc_samples
        and is_study_sample(meta.at[s, "sample_type"])
        and not dilution[s].strip()
    ]
    if args.qc_samples and len(qc_samples) < 3:
        print(
            f"[pcpfm] Only {len(qc_samples)} samples match --qc-samples '{args.qc_samples}'; "
            "QC CV and D-ratio need at least 3 and are skipped."
        )

    alignment = {}
    run_info = record.data["runs"].get("samples", {})
    stem_map = {x["stem"]: x["sample"] for x in record.data["inputs"] if x.get("stem")}
    for stem, info in run_info.get("alignment", {}).items():
        alignment[stem_map.get(stem, stem)] = info

    pca_result = pca(X)
    smet, corr = sample_metrics(X, meta, alignment, pca_result)
    fmet = feature_metrics(ft, X, meta, qc_samples, study_samples)
    summary = summarize(
        ft,
        smet,
        fmet,
        qc_samples,
        pca_result,
        record.step("blank_filter").get("summary"),
    )
    summary["study_samples"] = study_samples

    qc_dir = record.subdir("qc")
    smet.to_csv(os.path.join(qc_dir, "sample_qc.tsv"), sep="\t")
    fmet.to_csv(os.path.join(qc_dir, "feature_qc.tsv"), sep="\t", index=False)
    record.add_artifact(
        "sample_qc",
        os.path.join(qc_dir, "sample_qc.tsv"),
        "Per-sample QC metrics and flags",
        "qc",
    )
    record.add_artifact(
        "feature_qc",
        os.path.join(qc_dir, "feature_qc.tsv"),
        "Per-feature QC metrics",
        "qc",
    )

    # the filtered table keeps all samples of the input table, including blanks if any were there
    smet_all = smet.reindex(all_cols).fillna({"flags": ""})
    filtered, fcols, flog = filter_table(ft, all_cols, fmet, smet_all, args)
    fpath = os.path.join(record.subdir("tables"), "feature_table_filtered.tsv")
    filtered.to_csv(fpath, sep="\t", index=False)
    record.add_artifact(
        "feature_table_filtered",
        fpath,
        "Feature table after QC filtering: " + "; ".join(flog[1:-1]),
        "qc",
        n_features=len(filtered),
        n_samples=len(fcols),
        sample_columns=list(fcols),
        source=args.table,
    )
    summary["filtering"] = flog

    # analysis table: normalization, batch correction, imputation, log2; entry point for `pcpfm analyze`
    X_an, wlog = wrangle.prepare_analysis_table(
        filtered[fcols].astype(float), sample_info(record, fcols), args
    )
    analysis = filtered.copy()
    analysis[fcols] = X_an.round(4)
    apath = os.path.join(record.subdir("tables"), "feature_table_analysis.tsv")
    analysis.to_csv(apath, sep="\t", index=False)
    record.add_artifact(
        "feature_table_analysis",
        apath,
        "Filtered feature table after " + "; ".join(wlog) + "; input for pcpfm analyze",
        "qc",
        n_features=len(analysis),
        n_samples=len(fcols),
        sample_columns=list(fcols),
        source="feature_table_filtered",
        log2=not args.no_log2,
    )
    summary["wrangling"] = wlog

    html = report.write_report(
        os.path.join(qc_dir, "qc_report.html"),
        record,
        ft,
        X,
        meta,
        smet,
        fmet,
        corr,
        pca_result,
        summary,
        qc_samples,
    )
    record.add_artifact("qc_report", html, "QC report (HTML)", "qc")
    record.save()

    print(
        f"[pcpfm] QC: {summary['n_features']} features x {summary['n_samples']} samples; "
        f"{len(summary['flagged_samples'])} samples flagged."
    )
    for line in flog:
        print("[pcpfm]   " + line)
    print(f"[pcpfm] Report: {html}")
    return summary
