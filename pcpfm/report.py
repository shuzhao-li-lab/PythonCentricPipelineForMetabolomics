"""Generate self-contained HTML QC reports with figures embedded as PNG. In the spirit of MultiQC"""

import base64
import html
import io
import json
import os
from collections.abc import Sequence
from typing import Any

import matplotlib
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from numpy.typing import ArrayLike, NDArray

from .project import ProjectRecord

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import __version__

PALETTE = [
    "#4269d0",
    "#efb118",
    "#ff725c",
    "#6cc5b0",
    "#3ca951",
    "#ff8ab7",
    "#a463f2",
    "#97bbf5",
    "#9c6b4e",
]


def _colors(types: list[str]) -> tuple[list[str], dict[str, str]]:
    """Assign palette colors to sample types in first-occurrence order.

    Parameters
    ----------
    types : list[str]
        Sample-type labels in display order.

    Returns
    -------
    tuple[list[str], dict[str, str]]
        Colors aligned to types and the type-to-color mapping.
    """
    levels = list(dict.fromkeys(types))
    cmap = {t: PALETTE[i % len(PALETTE)] for i, t in enumerate(levels)}
    return [cmap[t] for t in types], cmap


def _legend(ax: Axes, cmap: dict[str, str]) -> None:
    """Add a sample-type color legend to an axes.

    Parameters
    ----------
    ax : Axes
        Axes to update.
    cmap : dict[str, str]
        Sample-type labels mapped to colors.
    """
    for t, c in cmap.items():
        ax.scatter([], [], color=c, label=t)
    ax.legend(fontsize=8, frameon=False, loc="best")


def _png(fig: Figure) -> str:
    """Encode a figure as PNG and close it.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure to encode and close.

    Returns
    -------
    str
        Base64-encoded PNG image.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def fig_per_sample(
    smet: pd.DataFrame, column: str, ylabel: str, log: bool = False
) -> str:
    """Plot a sample metric in injection order.

    Parameters
    ----------
    smet : pandas.DataFrame
        Sample QC metrics indexed by sample ID.
    column : str
        Sample-metric column to plot.
    ylabel : str
        Label for the vertical axis.
    log : bool, optional
        Whether to use a logarithmic y axis.

    Returns
    -------
    str
        Base64-encoded PNG image.
    """
    s = smet.sort_values("injection_order")
    colors, cmap = _colors(list(s["sample_type"]))
    fig, ax = plt.subplots(figsize=(10, 3.2))
    ax.bar(range(len(s)), s[column], color=colors)
    if log:
        ax.set_yscale("log")
    ax.set_xticks(range(len(s)))
    ax.set_xticklabels(s.index, rotation=90, fontsize=6)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Samples in injection order")
    _legend(ax, cmap)
    return _png(fig)


def fig_intensity_boxplot(X: pd.DataFrame, smet: pd.DataFrame) -> str:
    """Plot distributions of positive log2 intensities by sample.

    Parameters
    ----------
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    smet : pandas.DataFrame
        Sample QC metrics indexed by sample ID.

    Returns
    -------
    str
        Base64-encoded PNG image.
    """
    order = smet.sort_values("injection_order").index
    colors, cmap = _colors(list(smet.loc[order, "sample_type"]))
    data = [np.log2(X[s][X[s] > 0]) for s in order]
    fig, ax = plt.subplots(figsize=(10, 3.2))
    bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.7)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.7)
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=90, fontsize=6)
    ax.set_ylabel("log2 intensity (detected)")
    _legend(ax, cmap)
    return _png(fig)


def fig_correlation(corr: pd.DataFrame, smet: pd.DataFrame) -> str:
    """Plot sample correlations ordered by sample type and injection order.

    Parameters
    ----------
    corr : pandas.DataFrame
        Sample-by-sample correlation matrix.
    smet : pandas.DataFrame
        Sample QC metrics indexed by sample ID.

    Returns
    -------
    str
        Base64-encoded PNG image.
    """
    order = smet.sort_values(["sample_type", "injection_order"]).index
    c = corr.loc[order, order]
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(c.to_numpy(), cmap="viridis", vmin=np.nanmin(c.to_numpy()), vmax=1)
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=5)
    ax.set_yticklabels(order, fontsize=5)
    fig.colorbar(im, ax=ax, label="Spearman r (log2 intensity)")
    return _png(fig)


def fig_pca(
    pca_result: tuple[pd.DataFrame, NDArray[np.float64], pd.Series, float, int],
    smet: pd.DataFrame,
) -> str:
    """Plot the first two principal components colored by sample type.

    Parameters
    ----------
    pca_result : tuple[pandas.DataFrame, NDArray[np.float64], pandas.Series, float, int]
        PCA scores, variance ratios, T-squared values, threshold and component
        count, as returned by qc.pca. A result is required for this plot.
    smet : pandas.DataFrame
        Sample QC metrics indexed by sample ID.

    Returns
    -------
    str
        Base64-encoded PNG image.
    """
    scores, var, t2, limit, k = pca_result
    colors, cmap = _colors(list(smet.loc[scores.index, "sample_type"]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.scatter(
        scores["PC1"], scores["PC2"], c=colors, s=30, edgecolor="white", linewidth=0.5
    )
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)")
    for s in t2.index[t2 > limit]:
        ax.annotate(s, (scores.at[s, "PC1"], scores.at[s, "PC2"]), fontsize=7)
    _legend(ax, cmap)
    ax = axes[1]
    ax.bar(range(len(t2)), t2, color=colors)
    ax.axhline(limit, color="#c0392b", linestyle="--", linewidth=1, label="95% limit")
    ax.set_ylabel(f"Hotelling T² ({k} PCs)")
    ax.set_xticks(range(len(t2)))
    ax.set_xticklabels(t2.index, rotation=90, fontsize=5)
    return _png(fig)


def fig_feature_map(ft: pd.DataFrame) -> str:
    """Plot feature mass versus retention time with blank flags highlighted.

    Parameters
    ----------
    ft : pandas.DataFrame
        Feature table containing sample intensities and available feature
        annotations.

    Returns
    -------
    str
        Base64-encoded PNG image.
    """
    fig, ax = plt.subplots(figsize=(8, 4.5))
    flag = (
        ft["blank_flag"].astype(bool)
        if "blank_flag" in ft.columns
        else np.zeros(len(ft), bool)
    )
    size = (
        np.clip(np.log10(ft["peak_area"].clip(lower=1)) - 3, 0.5, 8)
        if "peak_area" in ft.columns
        else 3
    )
    ax.scatter(
        ft.loc[~flag, "rtime"],
        ft.loc[~flag, "mz"],
        s=size[~flag] if np.ndim(size) else size,
        color="#4269d0",
        alpha=0.5,
        label="feature",
        linewidth=0,
    )
    if flag.any():
        ax.scatter(
            ft.loc[flag, "rtime"],
            ft.loc[flag, "mz"],
            s=size[flag] if np.ndim(size) else size,
            color="#ff725c",
            alpha=0.7,
            label="blank-flagged",
            linewidth=0,
        )
    ax.set_xlabel("Retention time (s)")
    ax.set_ylabel("m/z")
    ax.legend(fontsize=8, frameon=False)
    return _png(fig)


def fig_hist(
    values: ArrayLike,
    xlabel: str,
    lines: Sequence[float] = (),
    bins: int = 50,
    xlim: tuple[float, float] | None = None,
) -> str:
    """Plot a histogram of finite values with optional reference lines.

    Parameters
    ----------
    values : ArrayLike
        Numeric observations; non-finite entries are omitted.
    xlabel : str
        Label for the horizontal axis.
    lines : Sequence[float], optional
        Positions of vertical reference lines.
    bins : int, optional
        Number of histogram bins.
    xlim : tuple[float, float] | None, optional
        Optional lower and upper x-axis bounds.

    Returns
    -------
    str
        Base64-encoded PNG image.
    """
    fig, ax = plt.subplots(figsize=(5, 3.2))
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if xlim:
        v = np.clip(v, *xlim)
    ax.hist(v, bins=bins, color="#4269d0", alpha=0.85)
    for x in lines:
        ax.axvline(x, color="#c0392b", linestyle="--", linewidth=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Features")
    return _png(fig)


def _table_html(df: pd.DataFrame, max_rows: int = 500) -> str:
    """Render a bounded preview of a DataFrame as an HTML table.

    Parameters
    ----------
    df : pandas.DataFrame
        Table to display.
    max_rows : int, optional
        Maximum number of rows to render.

    Returns
    -------
    str
        Escaped HTML table preview.
    """
    return df.head(max_rows).to_html(
        classes="data", border=0, na_rep="", float_format=lambda x: f"{x:.4g}"
    )


CSS = """
:root { --fg:#1d1d1f; --bg:#ffffff; --muted:#6e6e73; --line:#e5e5ea; --warn:#b3261e; --card:#f7f7f9; }
@media (prefers-color-scheme: dark) { :root { --fg:#f2f2f7; --bg:#1c1c1e; --muted:#a1a1a6; --line:#3a3a3c; --warn:#ff6b60; --card:#2c2c2e; } }
body { font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; color: var(--fg); background: var(--bg);
       max-width: 1150px; margin: 0 auto; padding: 16px; line-height: 1.45; }
h1 { font-size: 1.6em; margin-bottom: 0.2em; } h2 { margin-top: 1.8em; border-bottom: 1px solid var(--line); padding-bottom: 4px; }
.muted { color: var(--muted); font-size: 0.9em; }
.tiles { display: grid; grid-template-columns: repeat(auto-fill, minmax(170px, 1fr)); gap: 10px; }
.tile { background: var(--card); border-radius: 8px; padding: 10px 12px; }
.tile .v { font-size: 1.5em; font-weight: 600; } .tile .k { color: var(--muted); font-size: 0.85em; }
img { max-width: 100%; background: #fff; border-radius: 4px; }
.scroll { overflow-x: auto; max-height: 520px; }
table.data { border-collapse: collapse; font-size: 0.8em; }
table.data th, table.data td { border-bottom: 1px solid var(--line); padding: 3px 8px; text-align: right; white-space: nowrap; }
table.prov td { text-align: left; white-space: normal; vertical-align: top; }
table.prov td:first-child { white-space: nowrap; font-weight: 600; }
table.data th { position: sticky; top: 0; background: var(--bg); }
.warn { color: var(--warn); }
pre { background: var(--card); padding: 10px; border-radius: 6px; overflow-x: auto; font-size: 0.85em; }
"""


def write_report(
    path: str,
    record: ProjectRecord,
    ft: pd.DataFrame,
    X: pd.DataFrame,
    meta: pd.DataFrame,
    smet: pd.DataFrame,
    fmet: pd.DataFrame,
    corr: pd.DataFrame,
    pca_result: tuple[pd.DataFrame, NDArray[np.float64], pd.Series, float, int] | None,
    summary: dict[str, Any],
    qc_samples: list[str],
) -> str:
    """Write a self-contained HTML QC report with embedded PNG figures.

    Parameters
    ----------
    path : str
        File path to read or write.
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.
    ft : pandas.DataFrame
        Feature table containing sample intensities and available feature
        annotations.
    X : pandas.DataFrame
        Intensity-only DataFrame with features as rows and sample IDs as
        columns.
    meta : pandas.DataFrame
        Sample metadata indexed by sample ID.
    smet : pandas.DataFrame
        Sample QC metrics indexed by sample ID.
    fmet : pandas.DataFrame
        Feature QC metrics in feature-table row order.
    corr : pandas.DataFrame
        Sample-by-sample correlation matrix.
    pca_result : tuple[pandas.DataFrame, NDArray[np.float64], pandas.Series, float, int] | None
        PCA scores, variance ratios, T-squared values, threshold and component
        count, or None when unavailable.
    summary : dict[str, Any]
        Summary statistics and processing information.
    qc_samples : list[str]
        Sample IDs selected as technical QC replicates.

    Returns
    -------
    str
        The written report path.
    """
    e = html.escape
    tiles = [
        ("Features", summary["n_features"]),
        ("Samples", summary["n_samples"]),
        ("Median features / sample", int(summary["median_features_per_sample"])),
        ("Median detection rate", f"{summary['median_detection_rate']:.2f}"),
    ]
    if "blank_flagged_features" in summary:
        tiles.append(("Blank-flagged features", summary["blank_flagged_features"]))
    if summary.get("qc_cv_median") is not None:
        tiles.append(("Median QC CV %", summary["qc_cv_median"]))
        tiles.append(("Features QC CV < 30%", f"{summary['qc_cv_below_30_pct']}%"))
    if summary.get("dilution_r_above_0.8_pct") is not None:
        tiles.append(("Dilution r > 0.8", f"{summary['dilution_r_above_0.8_pct']}%"))
    tiles.append(("Flagged samples", len(summary["flagged_samples"])))

    parts = [
        (
            f"<h1>PCPFM QC report</h1><p class='muted'>Project {e(record.data['project']['name'])} · "
            f"pcpfm {__version__} · {e(record.data['project']['updated'])}</p>"
        ),
        "<div class='tiles'>"
        + "".join(
            f"<div class='tile'><div class='v'>{e(str(v))}</div><div class='k'>{e(k)}</div></div>"
            for k, v in tiles
        )
        + "</div>",
    ]
    if summary["flagged_samples"]:
        parts.append(
            "<h2>Flagged samples</h2><ul>"
            + "".join(
                f"<li><b>{e(s)}</b>: <span class='warn'>{e(f)}</span></li>"
                for s, f in summary["flagged_samples"].items()
            )
            + "</ul><p class='muted'>Flags are advisory: robust z-score beyond ±3 (median/MAD), "
            "Hotelling T² above the 95% limit, or failed Asari RT alignment. "
            "Designed variation, such as a dilution series, can trigger them. "
            "Flags compare each sample with others of the same sample_type and group; "
            "groups with fewer than 3 samples are not scored.</p>"
        )

    parts.append("<h2>Samples</h2>")
    parts.append(
        "<h3>Detected features</h3><img src='data:image/png;base64,"
        + fig_per_sample(smet, "n_features", "Detected features")
        + "'>"
    )
    parts.append(
        "<h3>Intensity distributions</h3><img src='data:image/png;base64,"
        + fig_intensity_boxplot(X, smet)
        + "'>"
    )
    parts.append(
        "<h3>Total intensity</h3><img src='data:image/png;base64,"
        + fig_per_sample(smet, "total_intensity", "Sum of intensities", log=True)
        + "'>"
    )
    if (
        "asari_landmark_peaks" in smet.columns
        and smet["asari_landmark_peaks"].notna().any()
    ):
        parts.append(
            "<h3>Asari RT alignment landmarks</h3><img src='data:image/png;base64,"
            + fig_per_sample(
                smet.fillna({"asari_landmark_peaks": 0}),
                "asari_landmark_peaks",
                "Landmark peaks",
            )
            + "'>"
            "<p class='muted'>Number of landmark peaks Asari used to align each sample to the reference; "
            "low numbers mean less reliable retention time alignment.</p>"
        )
    parts.append(
        "<h3>Sample correlation</h3><img src='data:image/png;base64,"
        + fig_correlation(corr, smet)
        + "'>"
    )
    if pca_result is not None:
        parts.append(
            "<h3>PCA</h3><img src='data:image/png;base64,"
            + fig_pca(pca_result, smet)
            + "'>"
            "<p class='muted'>Autoscaled log2 intensities of features detected in at least half of "
            "the samples; missing values set to half the feature minimum.</p>"
        )
    parts.append(
        "<h3>Sample metrics</h3><div class='scroll'>" + _table_html(smet) + "</div>"
    )

    parts.append("<h2>Features</h2>")
    parts.append(
        "<h3>Feature map</h3><img src='data:image/png;base64,"
        + fig_feature_map(ft)
        + "'>"
    )
    hists = [fig_hist(fmet["detection_rate"], "Detection rate", bins=40)]
    if "snr" in fmet.columns:
        hists.append(fig_hist(np.log10(fmet["snr"].clip(lower=1)), "log10 SNR (Asari)"))
    if "goodness_fitting" in fmet.columns:
        hists.append(
            fig_hist(fmet["goodness_fitting"], "Peak shape, goodness of fit (Asari)")
        )
    if "qc_cv" in fmet.columns:
        hists.append(
            fig_hist(
                fmet["qc_cv"],
                f"CV % in QC samples (n={len(qc_samples)})",
                lines=(20, 30),
                xlim=(0, 150),
            )
        )
    if "d_ratio" in fmet.columns:
        hists.append(
            fig_hist(
                fmet["d_ratio"],
                f"D-ratio % (QC SD / study SD, n={len(summary.get('study_samples', []))} study samples)",
                lines=(50,),
                xlim=(0, 200),
            )
        )
    if "dilution_r" in fmet.columns:
        hists.append(
            fig_hist(
                fmet["dilution_r"],
                "Pearson r, log intensity vs log concentration",
                lines=(0.8,),
            )
        )
    parts.append(
        "<div class='tiles' style='grid-template-columns: repeat(auto-fill, minmax(330px, 1fr))'>"
        + "".join(f"<img src='data:image/png;base64,{h}'>" for h in hists)
        + "</div>"
    )

    parts.append(
        "<h2>Filtering and wrangling</h2><pre>"
        + e("\n".join(summary.get("filtering", []) + summary.get("wrangling", [])))
        + "</pre>"
    )
    parts.append(
        "<h2>Provenance</h2><pre>"
        + e(
            json.dumps(
                {
                    "runs": {
                        k: {
                            x: v.get(x)
                            for x in (
                                "asari_version",
                                "n_files",
                                "output_dir",
                                "command",
                                "dropped_samples",
                            )
                        }
                        for k, v in record.data["runs"].items()
                    },
                    "steps": {
                        k: {x: v.get(x) for x in ("status", "finished", "summary")}
                        for k, v in record.data["steps"].items()
                    },
                },
                indent=2,
                default=str,
            )
        )
        + "</pre>"
    )

    doc = (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>PCPFM QC report</title><style>{CSS}</style></head><body>"
        + "\n".join(parts)
        + "</body></html>"
    )
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(doc)
    return path
