"""
``pcpfm replicates``: collapse technical replicates to one column per sample.

Run after ``pcpfm process`` when samples were injected in replicate.

1. Group sample columns by the metadata column given with --rep-col; samples with an
   empty value are not replicates and pass through unchanged (e.g. pooled QC samples).
2. Number replicates within a group by injection order (metadata row order without an
   injection_order column) and compute the Pearson correlation of log2(x + 1)
   intensities between each pair of replicates.
3. Flag groups whose minimum or maximum correlation is below the cutoff.
   max < cutoff: no replicate agrees with another and the sample is excluded.
   min < cutoff: replicates without any pair at or above the cutoff are dropped.
   --no-drop-flagged keeps all replicates of all groups.
4. Replace the remaining replicates by their per-feature median, stored in the column
   of the first remaining replicate (the representative); other replicate columns are
   dropped.

Outputs:

    replicates/replicate_correlations.xlsx   one row per replicate group
    tables/feature_table_replicates.tsv      collapsed table, raw intensities
"""

import itertools
import os
from argparse import Namespace
from typing import Any

import numpy as np
import pandas as pd

from .metadata import load_metadata
from .project import COMPLETED, FAILED, ProjectRecord
from .tables import load_table

STATUS_OK, STATUS_DROPPED, STATUS_EXCLUDED, STATUS_KEPT, STATUS_SINGLE = (
    "collapsed",
    "collapsed without flagged replicates",
    "excluded",
    "collapsed, flags ignored (--no-drop-flagged)",
    "single replicate",
)


def replicate_metadata(
    record: ProjectRecord, metadata_path: str | None, rep_col: str
) -> pd.DataFrame:
    """Load sample metadata holding the replicate assignment column.

    Parameters
    ----------
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.
    metadata_path : str | None
        Metadata CSV/TSV given on the command line, or None for the project copy.
    rep_col : str
        Metadata column with the replicate group of each sample.

    Returns
    -------
    pd.DataFrame
        String-valued metadata indexed by sample ID, with the original row order in
        ``_order``.

    Raises
    ------
    FileNotFoundError
        If no metadata is given and the project has no metadata copy.
    ValueError
        If the metadata lacks the replicate column.
    """
    if metadata_path:
        meta = load_metadata(metadata_path)
        source = metadata_path
    else:
        source = record.artifact_path("metadata", must_exist=False)
        if not source or not os.path.exists(source):
            raise FileNotFoundError(
                "No sample metadata: the project has no metadata copy (process was run "
                "without --metadata) and none was given with --metadata. Pass the "
                f"metadata file with the '{rep_col}' column via --metadata."
            )
        meta = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    if rep_col not in meta.columns:
        raise ValueError(
            f"Replicate column '{rep_col}' is not in the metadata {source}; "
            f"found {list(meta.columns)}"
        )
    meta["_order"] = range(len(meta))
    return meta.set_index("sample")


def replicate_groups(
    meta: pd.DataFrame, sample_cols: list[str], rep_col: str
) -> dict[str, list[str]]:
    """Group feature-table samples by their replicate assignment.

    Parameters
    ----------
    meta : pandas.DataFrame
        Metadata indexed by sample ID with the replicate column and ``_order``.
    sample_cols : list[str]
        Sample columns of the feature table.
    rep_col : str
        Metadata column with the replicate group of each sample.

    Returns
    -------
    dict[str, list[str]]
        Replicate samples by group in order of first appearance in the table. Samples
        within a group are ordered by injection_order, then metadata row order.

    Notes
    -----
    Samples absent from the metadata or with an empty group value are left out.
    Nonnumeric injection orders sort last.
    """
    meta = meta.reindex(sample_cols)
    group = meta[rep_col].fillna("").astype(str).str.strip()
    if "injection_order" in meta.columns:
        injection = pd.to_numeric(meta["injection_order"], errors="coerce")
    else:
        injection = pd.Series(np.nan, index=meta.index)
    order = pd.DataFrame(
        {
            "injection": injection,
            "row": meta["_order"].astype(float),
            "position": range(len(sample_cols)),
        },
        index=sample_cols,
    ).sort_values(["injection", "row", "position"], na_position="last")
    # groups in table order, members in injection order
    groups: dict[str, list[str]] = {g: [] for g in group[group != ""]}
    for s in order.index:
        if group[s]:
            groups[group[s]].append(s)
    return groups


def pair_label(i: int, j: int) -> str:
    """Name a replicate pair by 1-based replicate positions.

    Parameters
    ----------
    i, j : int
        0-based replicate positions within a group.

    Returns
    -------
    str
        Label such as "Rep1 vs Rep2".
    """
    return f"Rep{i + 1} vs Rep{j + 1}"


def assess_group(
    L: pd.DataFrame, members: list[str], cutoff: float, drop_flagged: bool
) -> dict[str, Any]:
    """Correlate the replicates of one group and decide which to keep.

    Parameters
    ----------
    L : pandas.DataFrame
        Intensities used for correlation, features as rows and samples as columns.
    members : list[str]
        Replicate samples of the group in replicate order.
    cutoff : float
        Minimum Pearson correlation for two replicates to agree.
    drop_flagged : bool
        Whether to drop replicates without any agreeing pair and exclude groups in
        which no pair agrees.

    Returns
    -------
    dict[str, Any]
        Pairwise correlations by pair label, min, max and mean correlation, flag text,
        kept and dropped replicates, representative sample and status.

    Notes
    -----
    A correlation is undefined (NaN) when a replicate has constant intensities, e.g. a
    failed injection with no detected features; it counts as below the cutoff.
    """
    pairs = {}
    agrees = {s: False for s in members}
    for (i, a), (j, b) in itertools.combinations(enumerate(members), 2):
        with np.errstate(invalid="ignore", divide="ignore"):
            r = float(np.corrcoef(L[a].to_numpy(), L[b].to_numpy())[0, 1])
        pairs[pair_label(i, j)] = r
        if r >= cutoff:
            agrees[a] = agrees[b] = True

    values = np.array(list(pairs.values()), dtype=float)
    finite = values[np.isfinite(values)]
    out = {
        "pairs": pairs,
        "min": float(finite.min()) if finite.size else np.nan,
        "max": float(finite.max()) if finite.size else np.nan,
        "mean": float(finite.mean()) if finite.size else np.nan,
        "flag": "",
    }
    if len(members) < 2:
        status, kept = STATUS_SINGLE, list(members)
    elif not any(agrees.values()):
        out["flag"] = f"max < {cutoff:g}"
        status, kept = STATUS_EXCLUDED, []
    elif not (values >= cutoff).all():
        out["flag"] = f"min < {cutoff:g}"
        status, kept = STATUS_DROPPED, [s for s in members if agrees[s]]
    else:
        status, kept = STATUS_OK, list(members)
    if out["flag"] and not drop_flagged:
        status, kept = STATUS_KEPT, list(members)
    out.update(
        status=status,
        kept=kept,
        dropped=[s for s in members if s not in kept],
        representative=kept[0] if kept else "",
    )
    return out


def collapse(
    ft: pd.DataFrame,
    sample_cols: list[str],
    groups: dict[str, list[str]],
    results: dict[str, dict[str, Any]],
) -> tuple[pd.DataFrame, list[str]]:
    """Replace each replicate group by the median of its kept replicates.

    Parameters
    ----------
    ft : pandas.DataFrame
        Feature table with annotation and sample columns.
    sample_cols : list[str]
        Sample columns of the feature table.
    groups : dict[str, list[str]]
        Replicate samples by group.
    results : dict[str, dict[str, Any]]
        Group assessments from assess_group.

    Returns
    -------
    tuple[pd.DataFrame, list[str]]
        Collapsed table and its sample columns in original table order. Missing
        intensities count as zero.
    """
    group_of = {s: g for g, members in groups.items() for s in members}
    X = ft[sample_cols].astype(float).fillna(0)
    out_cols = []
    for s in sample_cols:
        g = group_of.get(s)
        if g is None:
            out_cols.append(s)
        elif s == results[g]["representative"]:
            X[s] = X[results[g]["kept"]].median(axis=1)
            out_cols.append(s)
    lead = [c for c in ft.columns if c not in sample_cols]
    collapsed = pd.concat([ft[lead], X[out_cols]], axis=1)
    return collapsed, out_cols


def correlation_table(
    groups: dict[str, list[str]], results: dict[str, dict[str, Any]], rep_col: str
) -> pd.DataFrame:
    """Build the per-group replicate correlation report.

    Parameters
    ----------
    groups : dict[str, list[str]]
        Replicate samples by group, in replicate order.
    results : dict[str, dict[str, Any]]
        Group assessments from assess_group.
    rep_col : str
        Metadata column with the replicate group, used as the first column name.

    Returns
    -------
    pd.DataFrame
        One row per group: replicate sample names, pairwise correlations, min, max
        and average correlation, flag, status, dropped replicates and representative.
    """
    n_max = max((len(m) for m in groups.values()), default=0)
    rep_cols = [f"Rep{i + 1}" for i in range(n_max)]
    pair_cols = [pair_label(i, j) for i, j in itertools.combinations(range(n_max), 2)]
    rows = []
    for g, members in groups.items():
        res = results[g]
        row = {rep_col: g, "n_replicates": len(members)}
        row.update(dict(zip(rep_cols, members)))
        row.update(res["pairs"])
        row.update(
            {
                "min correlation": res["min"],
                "max correlation": res["max"],
                "average correlation": res["mean"],
                "flag": res["flag"],
                "status": res["status"],
                "dropped replicates": ", ".join(res["dropped"]),
                "representative": res["representative"],
            }
        )
        rows.append(row)
    columns = (
        [rep_col, "n_replicates"]
        + rep_cols
        + pair_cols
        + [
            "min correlation",
            "max correlation",
            "average correlation",
            "flag",
            "status",
            "dropped replicates",
            "representative",
        ]
    )
    df = pd.DataFrame(rows, columns=columns)
    numeric = pair_cols + ["min correlation", "max correlation", "average correlation"]
    df[numeric] = df[numeric].astype(float).round(4)
    return df


def run(args: Namespace, record: ProjectRecord) -> None:
    """Run the replicates step and record its status and summary.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for the ``replicates`` command.
    record : ProjectRecord
        Project record to update.

    Raises
    ------
    Exception
        Re-raises errors after printing them and marking the step as failed.
    """
    params = {k: v for k, v in vars(args).items() if k not in ("func", "argv")}
    record.start_step("replicates", params)
    try:
        summary = _replicates(args, record)
    except Exception as e:
        print(f"[pcpfm] Error: {e}")
        record.finish_step("replicates", FAILED, message=str(e))
        raise
    record.finish_step("replicates", COMPLETED, summary=summary)


def _replicates(args: Namespace, record: ProjectRecord) -> dict[str, Any]:
    """Assess replicate correlation, collapse replicates and register outputs.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for the ``replicates`` command.
    record : ProjectRecord
        Project record containing inputs, artifacts and pipeline state.

    Returns
    -------
    dict[str, Any]
        Step summary for the project record.

    Raises
    ------
    ValueError
        If the cutoff is outside [-1, 1] or no sample has a replicate group.
    """
    if not -1 <= args.pearson <= 1:
        raise ValueError(f"--pearson must be between -1 and 1, got {args.pearson}.")
    meta = replicate_metadata(record, args.metadata, args.rep_col)
    ft, sample_cols = load_table(record, args.table)
    groups = replicate_groups(meta, sample_cols, args.rep_col)
    if not groups:
        raise ValueError(
            f"No sample of table '{args.table}' has a value in column '{args.rep_col}'."
        )
    missing = [s for s in sample_cols if s not in meta.index]
    if missing:
        print(
            f"[pcpfm] Warning: {len(missing)} samples are not in the metadata and are "
            f"kept as they are: {missing[:10]}"
        )
    grouped = {s for members in groups.values() for s in members}
    ungrouped = [s for s in sample_cols if s not in grouped]
    if ungrouped:
        print(
            f"[pcpfm] {len(ungrouped)} samples have no replicate group and are kept as they are."
        )

    X = ft[sample_cols].astype(float).fillna(0)
    L = X if args.no_log2 else np.log2(X + 1)
    drop_flagged = not args.no_drop_flagged
    results = {
        g: assess_group(L, members, args.pearson, drop_flagged)
        for g, members in groups.items()
    }
    collapsed, out_cols = collapse(ft, sample_cols, groups, results)

    table = correlation_table(groups, results, args.rep_col)
    report_dir = record.subdir("replicates")
    report_path = os.path.join(report_dir, "replicate_correlations.xlsx")
    table.to_excel(report_path, sheet_name="replicate_correlations", index=False)
    record.add_artifact(
        "replicate_correlations",
        report_path,
        "Pairwise Pearson correlation between replicates, per replicate group",
        "replicates",
    )

    ft_path = os.path.join(record.subdir("tables"), "feature_table_replicates.tsv")
    collapsed.to_csv(ft_path, sep="\t", index=False)
    record.add_artifact(
        "feature_table_replicates",
        ft_path,
        f"Table '{args.table}' with replicates replaced by their median in the column "
        "of the representative replicate",
        "replicates",
        n_features=len(collapsed),
        n_samples=len(out_cols),
        sample_columns=list(out_cols),
        source_table=args.table,
        replicates={
            r["representative"]: r["kept"] for r in results.values() if r["kept"]
        },
    )
    record.save()

    flagged = [g for g, r in results.items() if r["flag"]]
    excluded = [g for g, r in results.items() if r["status"] == STATUS_EXCLUDED]
    n_dropped = sum(
        len(r["dropped"]) for r in results.values() if r["status"] == STATUS_DROPPED
    )
    print(
        f"[pcpfm] Replicates: {len(groups)} groups; {len(flagged)} flagged, "
        f"{len(excluded)} excluded, {n_dropped} replicates dropped; "
        f"{len(sample_cols)} -> {len(out_cols)} samples. Report: {record.relpath(report_path)}"
    )
    return {
        "n_groups": len(groups),
        "n_flagged": len(flagged),
        "n_excluded": len(excluded),
        "n_replicates_dropped": n_dropped,
        "flagged": flagged,
        "excluded": excluded,
        "n_samples_in": len(sample_cols),
        "n_samples_out": len(out_cols),
    }
