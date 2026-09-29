"""
``pcpfm analyze``: per-feature statistics against metadata variables, on study samples.

Input is the analysis table from ``pcpfm qc`` (log2 transformed by default).
Study samples exclude blanks, QC and standards by sample_type unless --samples is given.

For each variable (metadata column):
    categorical, 2 groups    Welch t-test; log2 fold change (difference of means of log2 values)
    categorical, >2 groups   one-way ANOVA; eta squared; largest difference between group means
    continuous               linear regression of feature intensity on the variable
P-values are adjusted by Benjamini-Hochberg (q). Figures: volcano plot, p-value histogram,
heatmap of top features, boxplots (or scatter plots) of top features.

Outputs go to analysis/<variable>/ and analysis/analysis_report.html.
"""

import base64
import html
import os
import re
import warnings

import matplotlib
import matplotlib.patches

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from . import __version__  # noqa: E402
from .project import COMPLETED, FAILED, now  # noqa: E402
from .qc import load_table, parse_selector, sample_info, select  # noqa: E402
from .report import CSS, PALETTE  # noqa: E402

# sample_type substrings that mark non-study samples
NON_STUDY_TYPES = ("blank", "qc", "pooled", "qstd", "std", "standard", "nist", "dda", "solvent", "process")
# metadata columns that are not study variables unless named explicitly
RESERVED = {"sample", "original_file_name", "file_path", "file", "run", "sample_type",
            "injection_order", "batch", "unspecified"}


# ---------------- samples and variables ----------------

def study_samples(meta, selector_text=None, exclude=()):
    if selector_text:
        samples = select(meta, parse_selector(selector_text))
    else:
        samples = [s for s in meta.index
                   if not any(t in str(meta.at[s, "sample_type"]).lower() for t in NON_STUDY_TYPES)]
    return [s for s in samples if s not in set(exclude)]


def natural_key(x):
    """Sort key that orders embedded numbers numerically: 'Veg(4:1)' < 'Veg(16:1)'."""
    return [(0, float(t)) if re.fullmatch(r"\d+(\.\d+)?", t) else (1, t.lower())
            for t in re.findall(r"\d+(?:\.\d+)?|[^\d]+", str(x))]


def _numeric(values):
    v = pd.to_numeric(values, errors="coerce")
    return v if v.notna().all() else None


def classify_variable(values, force=None, min_group_size=2):
    """
    values: Series over study samples (strings, "" = missing).
    Returns (kind, Series of usable values) with kind in {"categorical", "continuous", None}.
    """
    v = values[values.astype(str).str.strip() != ""].astype(str).str.strip()
    if len(v) < 4:
        return None, v
    num = _numeric(v)
    if force == "continuous" or (force is None and num is not None and num.nunique() >= 3):
        if num is None:
            num = pd.to_numeric(v, errors="coerce").dropna()
        return ("continuous", num) if num.nunique() >= 3 else (None, num)
    counts = v.value_counts()
    keep = counts[counts >= min_group_size].index
    v = v[v.isin(keep)]
    if len(keep) < 2 or (force is None and len(keep) > len(v) / 2):
        return None, v
    return "categorical", v


def select_variables(meta, samples, names=None, categorical=(), continuous=()):
    names = names or [c for c in meta.columns if c not in RESERVED and not c.startswith("_")]
    out = {}
    for c in names:
        if c not in meta.columns:
            raise ValueError(f"Variable '{c}' is not in the metadata.")
        force = "categorical" if c in categorical else "continuous" if c in continuous else None
        kind, values = classify_variable(meta.loc[samples, c].fillna(""), force)
        if kind:
            out[c] = (kind, values)
        else:
            print(f"[pcpfm] Skipping variable '{c}': not enough usable values or groups.")
    return out


# ---------------- statistics ----------------

def bh_fdr(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        res = np.empty(len(pv))
        res[order] = np.minimum(ranked, 1)
        q[ok] = res
    return q


def test_two_groups(Y, groups, reference=None):
    levels = sorted(groups.unique(), key=natural_key)
    if reference in levels:
        levels = [reference] + [x for x in levels if x != reference]
    a, b = (Y[groups.index[groups == lv]].to_numpy() for lv in levels)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t, p = stats.ttest_ind(b, a, axis=1, equal_var=False)
    res = pd.DataFrame(index=Y.index)
    res[f"mean_{levels[0]}"] = a.mean(axis=1)
    res[f"mean_{levels[1]}"] = b.mean(axis=1)
    res["log2FC"] = res[f"mean_{levels[1]}"] - res[f"mean_{levels[0]}"]
    res["t"], res["p"] = t, p
    return res, f"{levels[1]} vs {levels[0]}", "log2FC"


def test_anova(Y, groups):
    levels = sorted(groups.unique(), key=natural_key)
    arrays = [Y[groups.index[groups == lv]].to_numpy() for lv in levels]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")     # features constant within all groups give NaN
        F, p = stats.f_oneway(*arrays, axis=1)
    res = pd.DataFrame(index=Y.index)
    means = np.column_stack([a.mean(axis=1) for a in arrays])
    for lv, m in zip(levels, means.T):
        res[f"mean_{lv}"] = m
    grand = Y[groups.index].to_numpy().mean(axis=1, keepdims=True)
    ss_between = sum(a.shape[1] * (a.mean(axis=1) - grand[:, 0]) ** 2 for a in arrays)
    ss_total = ((Y[groups.index].to_numpy() - grand) ** 2).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        res["eta_squared"] = ss_between / ss_total
    res["max_diff"] = means.max(axis=1) - means.min(axis=1)
    res["F"], res["p"] = F, p
    return res, f"{len(levels)} groups", "max_diff"


def test_regression(Y, x):
    Yv = Y[x.index].to_numpy()
    xv = x.to_numpy(float)
    n = len(xv)
    xc = xv - xv.mean()
    yc = Yv - Yv.mean(axis=1, keepdims=True)
    sxx = (xc ** 2).sum()
    slope = (yc * xc).sum(axis=1) / sxx
    with np.errstate(divide="ignore", invalid="ignore"):
        r = (yc * xc).sum(axis=1) / np.sqrt(sxx * (yc ** 2).sum(axis=1))
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
    res = pd.DataFrame(index=Y.index)
    res["slope"] = slope
    res["intercept"] = Yv.mean(axis=1) - slope * xv.mean()
    res["r"] = r
    res["t"] = t
    res["p"] = 2 * stats.t.sf(np.abs(t), n - 2)
    return res, f"n={n}", "slope"


# ---------------- figures ----------------

def _save(fig, path):
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode()


def fig_volcano(res, effect, xlabel, alpha, path):
    fig, ax = plt.subplots(figsize=(6, 4.5))
    x, y = res[effect], -np.log10(res["p"].clip(lower=1e-300))
    sig = res["q"] < alpha
    ax.scatter(x[~sig], y[~sig], s=5, color="#9aa5b1", alpha=0.6, linewidth=0, label="not significant")
    ax.scatter(x[sig], y[sig], s=7, color="#ff725c", alpha=0.8, linewidth=0, label=f"q < {alpha}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("-log10 p")
    ax.legend(fontsize=8, frameon=False)
    return _save(fig, path)


def fig_pvalues(res, path):
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.hist(res["p"].dropna(), bins=40, color="#4269d0")
    ax.set_xlabel("p-value")
    ax.set_ylabel("Features")
    return _save(fig, path)


def _sample_order(values, kind):
    if kind == "continuous":
        return values.astype(float).sort_values().index
    return pd.Index(sorted(values.index, key=lambda i: (natural_key(values[i]), natural_key(i))))


def fig_heatmap(Y, top, values, kind, labels, path):
    order = _sample_order(values, kind)
    Z = Y.loc[top, order]
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1).replace(0, 1), axis=0)
    fig, (ax0, ax) = plt.subplots(2, 1, figsize=(10, 0.14 * len(top) + 2),
                                  gridspec_kw={"height_ratios": [0.4, max(len(top) / 10, 2)]}, sharex=True)
    if kind == "categorical":
        levels = list(dict.fromkeys(values[order]))
        colors = [PALETTE[i % len(PALETTE)] if len(levels) <= len(PALETTE) else plt.get_cmap("tab20")(i % 20)
                  for i in range(len(levels))]
        ax0.imshow([[levels.index(v) for v in values[order]]], aspect="auto",
                   cmap=matplotlib.colors.ListedColormap(colors), vmin=-0.5, vmax=len(levels) - 0.5)
        handles = [matplotlib.patches.Patch(color=c, label=lv) for lv, c in zip(levels, colors)]
        ax0.legend(handles=handles, loc="lower left", bbox_to_anchor=(1.01, 0), fontsize=6, frameon=False,
                   ncol=1 + len(levels) // 12)
    else:
        v = values[order].astype(float)
        ax0.imshow([v.to_numpy()], aspect="auto", cmap="viridis")
        ax0.set_title(f"samples ordered by value, {v.min():.3g} to {v.max():.3g}", fontsize=7, loc="left")
    ax0.set_yticks([])
    im = ax.imshow(Z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([labels[i] for i in top], fontsize=5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=5)
    fig.colorbar(im, ax=[ax0, ax], label="z-score", shrink=0.5)
    return _save(fig, path)


def fig_top_features(Y, top, values, kind, labels, path, ncol=10):
    n = len(top)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 2.0, nrow * 1.8), squeeze=False)
    if kind == "categorical":
        levels = sorted(values.unique(), key=natural_key)
        idx = {lv: values.index[values == lv] for lv in levels}
    x = values.astype(float) if kind == "continuous" else None
    for k, ax in enumerate(axes.flat):
        if k >= n:
            ax.axis("off")
            continue
        f = top[k]
        if kind == "categorical":
            data = [Y.loc[f, idx[lv]].to_numpy() for lv in levels]
            bp = ax.boxplot(data, showfliers=False, patch_artist=True, widths=0.6)
            for patch, c in zip(bp["boxes"], PALETTE * 5):
                patch.set_facecolor(c)
                patch.set_alpha(0.6)
            for i, d in enumerate(data):
                ax.scatter(np.full(len(d), i + 1) + np.random.default_rng(k).uniform(-0.15, 0.15, len(d)), d,
                           s=3, color="black", linewidth=0)
            ax.set_xticks(range(1, len(levels) + 1))
            ax.set_xticklabels(levels, rotation=90, fontsize=4)
        else:
            yv = Y.loc[f, x.index].to_numpy()
            ax.scatter(x, yv, s=4, color="#4269d0", linewidth=0)
            b, a = np.polyfit(x, yv, 1)
            xs = np.linspace(x.min(), x.max(), 10)
            ax.plot(xs, a + b * xs, color="#ff725c", linewidth=1)
            ax.tick_params(labelsize=4)
        ax.set_title(labels[f], fontsize=5)
        ax.tick_params(axis="y", labelsize=4)
    fig.tight_layout()
    return _save(fig, path)


# ---------------- entry point ----------------

def safe_name(name):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "variable"


def analyze_variable(Y, feat, name, kind, values, args, out_dir):
    if kind == "continuous" and name in (args.log_vars or []):
        if (values <= 0).any():
            raise ValueError(f"Cannot log-transform '{name}': non-positive values.")
        values = np.log2(values)
        name_label = f"log2({name})"
    else:
        name_label = name
    Yv = Y[values.index]
    Yv = Yv[Yv.std(axis=1) > 0]

    if kind == "continuous":
        res, desc, effect = test_regression(Yv, values)
        xlabel = f"slope (log2 intensity per unit {name_label})"
    elif values.nunique() == 2:
        res, desc, effect = test_two_groups(Yv, values, args.reference)
        xlabel = f"log2 fold change ({desc})"
    else:
        res, desc, effect = test_anova(Yv, values)
        xlabel = "largest difference between group means (log2)"
    res["q"] = bh_fdr(res["p"])
    res = feat.loc[res.index].join(res).sort_values("p")
    res = res.round({c: 4 for c in res.columns if c not in ("p", "q", "mz", "rtime")})

    os.makedirs(out_dir, exist_ok=True)
    res.to_csv(os.path.join(out_dir, "stats.tsv"), sep="\t", index=False)
    labels = {i: f"{feat.at[i, 'id_number']} {feat.at[i, 'mz']:.4f}@{feat.at[i, 'rtime']:.0f}s" for i in res.index}
    top = list(res.index[: args.top])
    figs = {
        "volcano": fig_volcano(res, effect, xlabel, args.alpha, os.path.join(out_dir, "volcano.png")),
        "pvalues": fig_pvalues(res, os.path.join(out_dir, "pvalue_histogram.png")),
        "heatmap": fig_heatmap(Yv, top[: args.heatmap_top], values, kind, labels, os.path.join(out_dir, "heatmap.png")),
        "top": fig_top_features(Yv, top, values, kind, labels, os.path.join(out_dir, "top_features.png")),
    }
    summary = {
        "variable": name,
        "label": name_label,
        "type": kind,
        "test": "linear regression" if kind == "continuous" else ("Welch t-test" if values.nunique() == 2 else "ANOVA"),
        "description": desc,
        "n_samples": int(len(values)),
        "levels": {str(k): int(v) for k, v in values.value_counts().items()} if kind == "categorical" else None,
        "n_features_tested": int(len(res)),
        "n_p_below_0.05": int((res["p"] < 0.05).sum()),
        f"n_q_below_{args.alpha}": int((res["q"] < args.alpha).sum()),
        "output_dir": out_dir,
    }
    return summary, figs, res


def run(args, record):
    params = {k: v for k, v in vars(args).items() if k != "func"}
    record.start_step("analyze", params)
    try:
        summary = _analyze(args, record)
    except Exception as e:
        record.finish_step("analyze", FAILED, message=str(e))
        raise
    record.finish_step("analyze", COMPLETED, summary=summary)


def _analyze(args, record):
    if not record.has_artifact(args.table):
        raise FileNotFoundError(f"Table '{args.table}' is not available; run `pcpfm qc` first.")
    if not record.data["artifacts"][args.table].get("log2", True):
        print("[pcpfm] Note: the table is not log2 transformed; statistics assume roughly normal intensities.")
    ft, sample_cols = load_table(record.artifact_path(args.table))
    meta = sample_info(record, sample_cols)
    exclude = []
    if args.exclude_flagged and record.has_artifact("sample_qc"):
        sq = pd.read_csv(record.artifact_path("sample_qc"), sep="\t", keep_default_na=False)
        exclude = sq.loc[sq["flags"] != "", "sample"].tolist()
    samples = study_samples(meta, args.samples, exclude)
    if len(samples) < 4:
        raise ValueError(f"Only {len(samples)} study samples found; use --samples column=value to select them.")
    print(f"[pcpfm] {len(samples)} study samples; excluded: "
          f"{sorted(set(sample_cols) - set(samples))[:20]}")

    variables = select_variables(meta, samples, args.variables, args.categorical or (), args.continuous or ())
    if not variables:
        raise ValueError("No usable metadata variables among study samples; name them with --variables.")

    Y = ft[sample_cols].astype(float)
    feat = ft[["id_number", "mz", "rtime"]]
    out_root = record.subdir("analysis")
    results, sections = [], []
    for name, (kind, values) in variables.items():
        print(f"[pcpfm] Analyzing '{name}' ({kind}, n={len(values)})")
        out_dir = os.path.join(out_root, safe_name(name))
        s, figs, res = analyze_variable(Y, feat, name, kind, values, args, out_dir)
        s["output_dir"] = record.relpath(out_dir)
        record.add_artifact(f"analysis_{safe_name(name)}", os.path.join(out_dir, "stats.tsv"),
                            f"{s['test']} of features vs {s['label']}", "analyze")
        results.append(s)
        sections.append((s, figs, res))
        print(f"[pcpfm]   {s['test']}: {s['n_features_tested']} features, "
              f"{s[f'n_q_below_{args.alpha}']} with q < {args.alpha}")

    table = pd.DataFrame([{k: v for k, v in s.items() if k != "levels"} for s in results])
    table.to_csv(os.path.join(out_root, "summary.tsv"), sep="\t", index=False)
    report_path = write_report(os.path.join(out_root, "analysis_report.html"), record, samples, sections, args)
    record.add_artifact("analysis_report", report_path, "Statistical analysis report (HTML)", "analyze")
    record.save()
    print(f"[pcpfm] Report: {report_path}")
    return {"n_study_samples": len(samples), "variables": results}


def write_report(path, record, samples, sections, args):
    e = html.escape
    parts = [f"<h1>PCPFM analysis report</h1><p class='muted'>Project {e(record.data['project']['name'])} · "
             f"pcpfm {__version__} · table {e(args.table)} · {len(samples)} study samples</p>"]
    parts.append("<p class='muted'><a href='#provenance'>Provenance</a></p>")
    parts.append("<table class='data'><tr><th>Variable</th><th>Type</th><th>Test</th><th>n</th>"
                 f"<th>Features</th><th>p &lt; 0.05</th><th>q &lt; {args.alpha}</th></tr>" + "".join(
                     f"<tr><td><a href='#{safe_name(s['variable'])}'>{e(s['label'])}</a></td><td>{s['type']}</td>"
                     f"<td>{e(s['test'])}</td><td>{s['n_samples']}</td><td>{s['n_features_tested']}</td>"
                     f"<td>{s['n_p_below_0.05']}</td><td>{s[f'n_q_below_{args.alpha}']}</td></tr>"
                     for s, _, _ in sections) + "</table>")
    for s, figs, res in sections:
        img = "<img src='data:image/png;base64,{}'>".format
        parts.append(f"<h2 id='{safe_name(s['variable'])}'>{e(s['label'])}</h2>"
                     f"<p class='muted'>{e(s['test'])}: {e(s['description'])}. Results: {e(s['output_dir'])}/stats.tsv</p>")
        if s["levels"]:
            parts.append("<p class='muted'>Groups: " + ", ".join(f"{e(k)} (n={v})" for k, v in s["levels"].items())
                         + "</p>")
        parts.append("<div class='tiles' style='grid-template-columns: repeat(auto-fill, minmax(420px, 1fr))'>"
                     + img(figs["volcano"]) + img(figs["pvalues"]) + "</div>")
        parts.append(f"<h3>Top {min(args.heatmap_top, len(res))} features</h3>" + img(figs["heatmap"]))
        parts.append(f"<h3>Top {min(args.top, len(res))} features</h3>" + img(figs["top"]))
        parts.append("<div class='scroll'>" + res.head(100).to_html(classes="data", border=0, index=False,
                                                                    float_format=lambda x: f"{x:.4g}") + "</div>")
    parts.append(provenance_html(record, samples, sections, args))
    doc = ("<!doctype html><html lang='en'><head><meta charset='utf-8'>"
           "<meta name='viewport' content='width=device-width, initial-scale=1'>"
           f"<title>PCPFM analysis report</title><style>{CSS}</style></head><body>" + "\n".join(parts)
           + "</body></html>")
    with open(path, "w", encoding="utf-8") as f:
        f.write(doc)
    return path


# ---------------- provenance ----------------

def _fn(*funcs):
    """Function names for display; strings are used as given (public names of library functions)."""
    names = [f if isinstance(f, str) else f"{f.__module__}.{f.__qualname__}" for f in funcs]
    return "<br>".join(f"<code>{html.escape(n)}</code>" for n in names)


def _command_text(argv_text):
    """Show the program as `pcpfm`, whether it was run as a script, module or entry point."""
    head, _, rest = argv_text.partition(" ")
    if head.endswith(("pcpfm", "__main__.py", "cli.py")):
        return "pcpfm " + rest
    return argv_text


def provenance_html(record, samples, sections, args):
    """Main steps from the project record, with the functions and parameters used."""
    import platform
    import sys

    import matplotlib as mpl
    import scipy

    from . import asari_runner, blanks, metadata, qc, wrangle

    e = html.escape
    d = record.data
    rows = []

    meta = d.get("metadata") or {}
    rows.append(("Sample metadata",
                 f"{meta.get('n_samples', 'no')} samples from {e(str(meta.get('source', 'none')))}; "
                 f"{sum(1 for x in d['inputs'] if x.get('file'))} data files resolved",
                 _fn(metadata.load_metadata, metadata.resolve_inputs)))
    for name, run in d["runs"].items():
        rows.append((f"Feature extraction, Asari run '{e(name)}'",
                     f"Asari {e(str(run.get('asari_version')))}, {run.get('n_files')} files; "
                     f"<code>{e(str(run.get('command', '')))}</code>; "
                     f"dropped samples: {e(str(run.get('dropped_samples') or 'none'))}",
                     _fn(asari_runner.run_asari)))
    proc = record.step("process").get("params", {})
    if proc:
        rows.append(("Feature table",
                     f"Asari {e(str(proc.get('asari_table')))} table of the sample run; "
                     "columns renamed from file names to sample IDs",
                     _fn(asari_runner.feature_table_path)))
    bf = record.step("blank_filter")
    if bf.get("status") == "completed":
        b = bf.get("summary", {})
        detail = (f"mode {e(str(b.get('mode')))}; flagged {b.get('n_flagged')} of {b.get('n_features')} features "
                  f"where blank mean x {b.get('ratio')} >= sample mean (non-zero means)")
        if b.get("mode") == "separate_runs":
            detail += (f"; matched by closest m/z and rtime within {b.get('ppm')} ppm and {b.get('rt_tol')} s "
                       f"after RT offset {b.get('rt_offset_seconds')} s ({b.get('rt_offset_pairs')} pairs); "
                       f"{b.get('n_matched')} features matched")
            funcs = _fn(blanks.estimate_rt_offset, blanks.match_features, blanks.flag_blank_features)
        else:
            funcs = _fn(blanks.blank_filter_same_run, blanks.flag_blank_features)
        rows.append(("Blank filtering", detail, funcs))
    else:
        rows.append(("Blank filtering", e(str(bf.get("status"))) + " " + e(str(bf.get("message", ""))), ""))

    q = record.step("qc")
    qs = q.get("summary", {})
    qp = q.get("params", {})
    rows.append(("Quality control",
                 f"QC samples selected by <code>{e(str(qp.get('qc_samples')))}</code>: "
                 f"{e(', '.join(qs.get('qc_samples', [])) or 'none')}; "
                 f"flagged samples: {e(', '.join(qs.get('flagged_samples', {})) or 'none')}",
                 _fn(qc.sample_metrics, qc.feature_metrics, qc.pca)))
    rows.append(("Filtering", "<br>".join(e(x) for x in qs.get("filtering", [])), _fn(qc.filter_table)))
    rows.append(("Normalization, batch correction, imputation, transform",
                 "<br>".join(e(x) for x in qs.get("wrangling", [])),
                 _fn(wrangle.normalize, wrangle.batch_correct, wrangle.impute, wrangle.log2)))

    rows.append(("Study samples",
                 f"{len(samples)} samples, selected by "
                 + (f"<code>{e(args.samples)}</code>" if args.samples else
                    "sample_type not matching " + e(", ".join(NON_STUDY_TYPES)))
                 + ("; samples flagged by QC excluded" if args.exclude_flagged else ""),
                 _fn(study_samples, select_variables)))
    test_funcs = {"Welch t-test": (test_two_groups, "scipy.stats.ttest_ind"),
                  "ANOVA": (test_anova, "scipy.stats.f_oneway"),
                  "linear regression": (test_regression, "scipy.stats.t.sf")}
    for s, _, _ in sections:
        extra = {"Welch t-test": "scipy.stats.ttest_ind(equal_var=False); log2FC = difference of group means",
                 "ANOVA": "scipy.stats.f_oneway; eta squared = SS between / SS total",
                 "linear regression": "ordinary least squares per feature; two-sided p from t distribution"}[s["test"]]
        rows.append((f"Statistics: {e(s['label'])}",
                     f"{e(s['test'])}"
                     + (f" ({e(s['description'])})" if s["type"] == "categorical" else "")
                     + f", n = {s['n_samples']}, "
                     f"{s['n_features_tested']} features with non-zero variance; {e(extra)}; "
                     f"Benjamini-Hochberg q-values, threshold {args.alpha}",
                     _fn(*test_funcs[s["test"]], bh_fdr)))
    rows.append(("Figures",
                 f"volcano plot, p-value histogram; heatmap of top {args.heatmap_top} features "
                 f"(z-scored per feature); boxplots or scatter plots of top {args.top} features by p-value",
                 _fn(fig_volcano, fig_pvalues, fig_heatmap, fig_top_features)))

    table = ("<table class='data prov'><tr><th>Step</th><th>What was done</th><th>Functions</th></tr>"
             + "".join(f"<tr><td>{a}</td><td>{b}</td><td>{c}</td></tr>" for a, b, c in rows) + "</table>")

    asari_versions = sorted({str(r.get("asari_version")) for r in d["runs"].values()})
    versions = (f"pcpfm {__version__}; Asari {', '.join(asari_versions) or 'n/a'}; Python {sys.version.split()[0]}; "
                f"numpy {np.__version__}; pandas {pd.__version__}; scipy {scipy.__version__}; "
                f"matplotlib {mpl.__version__}; {platform.platform()}")
    # latest completed command of each subcommand, i.e. the ones that produced the current outputs
    latest = {}
    for h in d["history"]:
        cmd = _command_text(h["command"])
        sub = cmd.split()[1] if len(cmd.split()) > 1 else cmd
        if h.get("status") == "completed" and sub != "analyze":
            latest[sub] = (h["time"], cmd)
    lines = [f"{t}  {c}" for t, c in sorted(latest.values())]
    lines.append(f"{now()}  {_command_text(' '.join(sys.argv))}   (this report)")
    commands = "\n".join(lines)
    return ("<h2 id='provenance'>Provenance</h2>"
            f"<p class='muted'>From the project record {e(record.path)}.</p>"
            "<div class='scroll'>" + table + "</div>"
            f"<h3>Software</h3><p>{e(versions)}</p>"
            f"<h3>Commands</h3><pre>{e(commands)}</pre>")
