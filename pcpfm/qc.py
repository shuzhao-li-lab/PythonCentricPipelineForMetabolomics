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

import os

import numpy as np
import pandas as pd
from scipy import stats

from .blanks import BLANK_COLUMNS
from .metadata import is_blank
from .process import FEATURE_COLUMNS
from .project import COMPLETED, FAILED

ZSCORE_LIMIT = 3.0


def robust_z(x):
    x = pd.Series(x, dtype=float)
    mad = 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))
    if not mad or np.isnan(mad):
        return x * 0.0
    return (x - np.nanmedian(x)) / mad


def rowwise_pearson(Y, x, min_n=4):
    """Pearson r of each row of Y (NaN = missing) with vector x, over pairwise complete values."""
    mask = ~np.isnan(Y)
    n = mask.sum(axis=1)
    Xb = np.where(mask, x, 0.0)
    Yz = np.where(mask, Y, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        xm = Xb.sum(axis=1) / n
        ym = Yz.sum(axis=1) / n
        dx = np.where(mask, Xb - xm[:, None], 0.0)
        dy = np.where(mask, Yz - ym[:, None], 0.0)
        r = (dx * dy).sum(axis=1) / np.sqrt((dx ** 2).sum(axis=1) * (dy ** 2).sum(axis=1))
    r[n < min_n] = np.nan
    return r


def parse_selector(text):
    """'sample_type=pooled' -> ('sample_type', {'pooled'}); values separated by commas."""
    if not text:
        return None
    key, _, values = text.partition("=")
    return key.strip(), {v.strip().lower() for v in values.split(",") if v.strip()}


# ---------------- loading ----------------

def load_table(path):
    ft = pd.read_csv(path, sep="\t")
    ft["id_number"] = ft["id_number"].astype(str)
    known = set(FEATURE_COLUMNS) | set(BLANK_COLUMNS)
    sample_cols = [c for c in ft.columns if c not in known]
    return ft, sample_cols


def sample_info(record, sample_cols):
    """Metadata rows for the table's samples, in table order, with injection_order filled."""
    path = record.artifact_path("metadata", must_exist=False)
    if path and os.path.exists(path):
        meta = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        meta["_order"] = range(len(meta))
        meta = meta.set_index("sample").reindex(sample_cols)
    else:
        meta = pd.DataFrame(index=pd.Index(sample_cols, name="sample"))
        meta["_order"] = range(len(meta))
    meta["sample_type"] = meta.get("sample_type", pd.Series("", index=meta.index)).fillna("")
    meta.loc[meta["sample_type"] == "", "sample_type"] = "unspecified"
    if "injection_order" in meta.columns:
        meta["injection_order"] = pd.to_numeric(meta["injection_order"], errors="coerce")
    else:
        meta["injection_order"] = meta["_order"]
    return meta.drop(columns="_order")


def select(meta, selector):
    if selector is None:
        return []
    key, values = selector
    if key not in meta.columns:
        return []
    return [s for s in meta.index if str(meta.at[s, key]).lower() in values]


# ---------------- computations ----------------

def pca(X, min_detection=0.5, max_components=10):
    """
    PCA on log2 intensities of features detected in >= min_detection of samples.
    Zeros are replaced by half the feature minimum; data are autoscaled.
    X: features x samples. Returns (scores DataFrame, explained variance ratio, T2 Series, T2 limit, k).
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
    var = S ** 2 / np.sum(S ** 2)
    scores = U[:, :ncomp] * S[:ncomp]
    # Hotelling T2 on the components explaining 80% of variance (2 to 5 components)
    k = int(np.clip(np.searchsorted(np.cumsum(var), 0.8) + 1, 2, min(5, ncomp)))
    lam = S[:k] ** 2 / (n - 1)
    t2 = np.sum(scores[:, :k] ** 2 / lam, axis=1)
    limit = k * (n - 1) / (n - k) * stats.f.ppf(0.95, k, n - k) if n > k else np.inf
    cols = [f"PC{i + 1}" for i in range(ncomp)]
    return (pd.DataFrame(scores, index=X.columns, columns=cols), var[:ncomp],
            pd.Series(t2, index=X.columns), float(limit), k)


def sample_metrics(X, meta, alignment=None, pca_result=None):
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
        scores, _, t2, limit, _ = pca_result
        df["PC1"], df["PC2"] = scores["PC1"].round(3), scores["PC2"].round(3)
        df["pca_T2"] = t2.round(2)
    if alignment:
        df["asari_landmark_peaks"] = [alignment.get(s, {}).get("landmark_peaks", np.nan) for s in df.index]
        df["asari_rt_aligned"] = [alignment.get(s, {}).get("aligned", True) for s in df.index]

    df["z_n_features"] = robust_z(df["n_features"]).round(2)
    df["z_total_intensity"] = robust_z(np.log2(df["total_intensity"] + 1)).round(2)
    df["z_median_correlation"] = robust_z(df["median_correlation"]).round(2)

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


def feature_metrics(ft, X, meta, qc_samples, study_samples):
    cols = [c for c in ("id_number", "mz", "rtime", "snr", "goodness_fitting", "cSelectivity", "blank_flag")
            if c in ft.columns]
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


def summarize(ft, smet, fmet, qc_samples, pca_result, blank_summary):
    s = {
        "n_features": int(len(ft)),
        "n_samples": int(len(smet)),
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
        s["qc_cv_below_20_pct"] = round(float((cv < 20).mean() * 100), 1) if len(cv) else None
        s["qc_cv_below_30_pct"] = round(float((cv < 30).mean() * 100), 1) if len(cv) else None
    if "d_ratio" in fmet.columns:
        d = fmet["d_ratio"].dropna()
        s["d_ratio_below_50_pct"] = round(float((d < 50).mean() * 100), 1) if len(d) else None
    if "dilution_r" in fmet.columns:
        r = fmet["dilution_r"].dropna()
        s["dilution_features_tested"] = int(len(r))
        s["dilution_r_above_0.8_pct"] = round(float((r > 0.8).mean() * 100), 1) if len(r) else None
    if pca_result is not None:
        s["pca_explained_variance"] = [round(float(v), 4) for v in pca_result[1][:5]]
        s["pca_T2_limit_95"] = round(pca_result[3], 2)
        s["pca_T2_components"] = pca_result[4]
    return s


# ---------------- wrangling ----------------

def filter_table(ft, sample_cols, fmet, smet, args):
    """Remove features and, optionally, flagged samples; intensities unchanged. Returns (table, sample_cols, log)."""
    log = [f"start: {len(ft)} features, {len(sample_cols)} samples"]
    keep = pd.Series(True, index=ft.index)
    if not args.keep_blank_features and "blank_flag" in ft.columns:
        keep &= ~ft["blank_flag"].astype(bool)
        log.append(f"remove blank-flagged features: {int((~keep).sum())} removed")
    before = int(keep.sum())
    keep &= fmet["detection_rate"].to_numpy() >= args.min_detection
    keep &= fmet["detection_rate"].to_numpy() > 0
    log.append(f"detection rate >= {args.min_detection}: {before - int(keep.sum())} removed")
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

def run(args, record):
    params = {k: v for k, v in vars(args).items() if k != "func"}
    record.start_step("qc", params)
    try:
        summary = _qc(args, record)
    except Exception as e:
        record.finish_step("qc", FAILED, message=str(e))
        raise
    record.finish_step("qc", COMPLETED, summary=summary)


def _qc(args, record):
    from . import report, wrangle

    ft, sample_cols = load_table(record.artifact_path(args.table))
    X = ft[sample_cols].astype(float).fillna(0)
    meta = sample_info(record, sample_cols)
    blank_like = [s for s in sample_cols if is_blank(meta.at[s, "sample_type"])]
    if blank_like:
        print(f"[pcpfm] Excluding {len(blank_like)} blank samples from QC metrics.")
        sample_cols = [s for s in sample_cols if s not in blank_like]
        X, meta = X[sample_cols], meta.loc[sample_cols]

    qc_samples = select(meta, parse_selector(args.qc_samples))
    study_samples = [s for s in sample_cols if s not in qc_samples]
    if args.qc_samples and len(qc_samples) < 3:
        print(f"[pcpfm] Only {len(qc_samples)} samples match --qc-samples '{args.qc_samples}'; "
              "QC CV and D-ratio need at least 3 and are skipped.")

    alignment = {}
    run_info = record.data["runs"].get("samples", {})
    stem_map = {x["stem"]: x["sample"] for x in record.data["inputs"] if x.get("stem")}
    for stem, info in run_info.get("alignment", {}).items():
        alignment[stem_map.get(stem, stem)] = info

    pca_result = pca(X)
    smet, corr = sample_metrics(X, meta, alignment, pca_result)
    fmet = feature_metrics(ft, X, meta, qc_samples, study_samples)
    summary = summarize(ft, smet, fmet, qc_samples, pca_result, record.step("blank_filter").get("summary"))

    qc_dir = record.subdir("qc")
    smet.to_csv(os.path.join(qc_dir, "sample_qc.tsv"), sep="\t")
    fmet.to_csv(os.path.join(qc_dir, "feature_qc.tsv"), sep="\t", index=False)
    record.add_artifact("sample_qc", os.path.join(qc_dir, "sample_qc.tsv"), "Per-sample QC metrics and flags", "qc")
    record.add_artifact("feature_qc", os.path.join(qc_dir, "feature_qc.tsv"), "Per-feature QC metrics", "qc")

    # the filtered table keeps all samples of the input table, including blanks if any were there
    all_cols = [c for c in load_table(record.artifact_path(args.table))[1]]
    smet_all = smet.reindex(all_cols).fillna({"flags": ""})
    filtered, fcols, flog = filter_table(ft, all_cols, fmet, smet_all, args)
    fpath = os.path.join(record.subdir("tables"), "feature_table_filtered.tsv")
    filtered.to_csv(fpath, sep="\t", index=False)
    record.add_artifact("feature_table_filtered", fpath, "Feature table after QC filtering: " + "; ".join(flog[1:-1]),
                        "qc", n_features=int(len(filtered)), n_samples=len(fcols), source=args.table)
    summary["filtering"] = flog

    # analysis table: normalization, batch correction, imputation, log2; entry point for `pcpfm analyze`
    X_an, wlog = wrangle.prepare_analysis_table(filtered[fcols].astype(float), sample_info(record, fcols), args)
    analysis = filtered.copy()
    analysis[fcols] = X_an.round(4)
    apath = os.path.join(record.subdir("tables"), "feature_table_analysis.tsv")
    analysis.to_csv(apath, sep="\t", index=False)
    record.add_artifact("feature_table_analysis", apath,
                        "Filtered feature table after " + "; ".join(wlog) + "; input for pcpfm analyze",
                        "qc", n_features=int(len(analysis)), n_samples=len(fcols), source="feature_table_filtered",
                        log2=not args.no_log2)
    summary["wrangling"] = wlog

    html = report.write_report(
        os.path.join(qc_dir, "qc_report.html"), record, ft, X, meta, smet, fmet, corr, pca_result, summary, qc_samples,
    )
    record.add_artifact("qc_report", html, "QC report (HTML)", "qc")
    record.save()

    print(f"[pcpfm] QC: {summary['n_features']} features x {summary['n_samples']} samples; "
          f"{len(summary['flagged_samples'])} samples flagged.")
    for line in flog:
        print("[pcpfm]   " + line)
    print(f"[pcpfm] Report: {html}")
    return summary
