"""
``pcpfm process``: mzML files -> feature tables.

1. Resolve input files from metadata and/or --files / --input.
2. Run Asari with --files: blanks and other samples in separate runs when metadata
   identifies blanks; otherwise a single run.
3. Rename sample columns from file names to sample IDs, and flag blank features.
4. Register everything in the project record.
"""

import os
from argparse import Namespace
from collections.abc import Sequence
from typing import Any

import pandas as pd

from . import asari_runner, blanks
from .metadata import list_input_files, load_metadata, resolve_inputs
from .project import COMPLETED, FAILED, SKIPPED, ProjectRecord
from .tables import load_asari_table


def _run_or_reuse(
    record: ProjectRecord,
    args: Namespace,
    cmd: list[str],
    run_name: str,
    files: Sequence[str],
    extra: Sequence[str],
    version: str,
) -> dict[str, Any]:
    """Run Asari unless an identical completed run is already recorded.

    A run is reusable when its Asari version, input files, parameters, and
    recorded output directory still match the requested run.

    Parameters
    ----------
    record : ProjectRecord
        Project record to read and update.
    args : argparse.Namespace
        Process command options, including mode, ppm, cores, and force.
    cmd : list of str
        Asari command and executable arguments.
    run_name : str
        Record key and output name for this run.
    files : sequence of str
        mzML files to process.
    extra : sequence of str
        Additional Asari arguments.
    version : str
        Detected Asari version.

    Returns
    -------
    dict[str, Any]
        The reused or newly recorded Asari run information.
    """
    sig = asari_runner.run_signature(
        files, list(extra) + [args.mode, str(args.ppm)], version
    )
    prev = record.data["runs"].get(run_name)
    usable = (
        prev
        and prev.get("status") == COMPLETED
        and os.path.exists(record.abspath(prev["output_dir"]))
    )
    if not args.force and usable and prev.get("signature") == sig:
        print(
            f"[pcpfm] Reusing Asari run '{run_name}' ({prev['output_dir']}); use --force to rerun."
        )
        return prev
    # tell the user why a run they may expect to be reused is redone
    if not args.force and usable:
        old = prev.get("asari_version")
        reason = (
            f"Asari version changed ({old} -> {version})"
            if old and old != version
            else "input files or parameters changed"
        )
        print(f"[pcpfm] Rerunning Asari '{run_name}': {reason}.")
    info = asari_runner.run_asari(
        cmd,
        run_name,
        files,
        record.subdir("asari"),
        args.mode,
        extra_args=extra,
        cores=args.cores,
        ppm=args.ppm,
        asari_version=version,
    )
    for k in ("files_list", "log", "output_dir"):
        info[k] = record.relpath(info[k])
    info["status"] = COMPLETED
    record.data["runs"][run_name] = info
    record.save()
    return info


def run(args: Namespace, record: ProjectRecord) -> None:
    """Run the process step and record its status and summary.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for the ``process`` command.
    record : ProjectRecord
        Project record to update.

    Raises
    ------
    Exception
        Re-raises processing errors after marking the step as failed.
    """
    params = {k: v for k, v in vars(args).items() if k not in ("func", "argv")}
    record.start_step("process", params)
    try:
        summary = _process(args, record)
    except Exception as e:
        record.finish_step("process", FAILED, message=str(e))
        raise
    record.finish_step("process", COMPLETED, summary=summary)
    for line in record.status_lines():
        print(line)


def _process(args: Namespace, record: ProjectRecord) -> dict[str, Any]:
    """Resolve inputs, run Asari, filter blanks, and register output tables.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for the ``process`` command.
    record : ProjectRecord
        Project record receiving metadata, run, and table artifacts.

    Returns
    -------
    dict[str, Any]
        A summary containing input, run, feature, and blank-filter counts.

    Raises
    ------
    ValueError
        If inputs are missing, unresolved, or otherwise unsuitable for
        processing.
    """
    # ---------- inputs ----------
    metadata = load_metadata(args.metadata) if args.metadata else None
    candidates = list_input_files(args.files, args.input)
    if metadata is None and not candidates:
        raise ValueError(
            "No input: give --metadata with a file_path column, or --files / --input."
        )
    if metadata is not None and "file_path" not in metadata.columns and not candidates:
        raise ValueError(
            "Metadata has no file_path column; also give --files or --input to locate mzML files."
        )

    inputs, unlisted = resolve_inputs(metadata, candidates, args.metadata)
    missing = [x["sample"] for x in inputs if not x["file"]]
    if missing:
        len_missing = len(missing)
        example = missing[:10] if len_missing > 10 else missing[:len_missing]
        print(
            f"[pcpfm] Warning: no data file found for {len_missing} samples: {example}"
        )
    if unlisted:
        len_unlisted = len(unlisted)
        example = unlisted[:10] if len_unlisted > 10 else unlisted[:len_unlisted]
        print(
            f"[pcpfm] Note: {len(unlisted)} mzML files are not in the metadata and are skipped: "
            f"{[os.path.basename(p) for p in example]}"
        )
    available = [x for x in inputs if x["file"]]
    if not available:
        raise ValueError("No data files to process.")

    separate = (
        not args.single_run
        and any(x["blank"] for x in available)
        and any(not x["blank"] for x in available)
    )
    for x in inputs:
        x["run"] = (
            None
            if not x["file"]
            else ("blanks" if separate and x["blank"] else "samples")
        )
    record.data["inputs"] = inputs

    if metadata is not None:
        meta = metadata.merge(
            pd.DataFrame(inputs)[["sample", "file", "run"]], on="sample", how="left"
        )
        meta_path = os.path.join(record.project_dir, "metadata.tsv")
        meta.to_csv(meta_path, sep="\t", index=False)
        record.data["metadata"] = {
            "source": os.path.abspath(args.metadata),
            "path": record.relpath(meta_path),
            "n_samples": len(metadata),
            "columns": list(metadata.columns),
        }
        record.add_artifact(
            "metadata",
            meta_path,
            "Sample metadata with resolved file paths and Asari run",
            "process",
        )
    record.save()

    # ---------- asari ----------
    cmd = asari_runner.asari_command(args.asari_cmd)
    version = asari_runner.check_asari(cmd)
    print(f"[pcpfm] Using Asari {version}: {' '.join(cmd)}")
    extra = args.asari_args.split() if args.asari_args else []
    stem_to_sample = {x["stem"]: x["sample"] for x in available}

    runs = {}
    for run_name in ["blanks", "samples"] if separate else ["samples"]:
        files = [x["file"] for x in available if x["run"] == run_name]
        runs[run_name] = _run_or_reuse(
            record, args, cmd, run_name, files, extra, version
        )
    for run_name in list(record.data["runs"]):
        if run_name not in runs:  # e.g. switching from two runs to --single-run
            del record.data["runs"][run_name]

    # ---------- tables ----------
    tables_dir = record.subdir("tables")
    sample_ft, sample_cols = load_asari_table(
        asari_runner.feature_table_path(
            record.abspath(runs["samples"]["output_dir"]), args.asari_table
        ),
        stem_to_sample,
    )
    blank_names = {x["sample"] for x in available if x["blank"]}

    record.start_step(
        "blank_filter",
        {"ratio": args.blank_ratio, "ppm": args.match_ppm, "rt_tol": args.match_rt},
    )
    if separate:
        # the blank table uses the full table: more sensitive detection gives more conservative filtering
        blank_ft, blank_cols = load_asari_table(
            asari_runner.feature_table_path(
                record.abspath(runs["blanks"]["output_dir"]), "full"
            ),
            stem_to_sample,
        )
        ft, bsum = blanks.blank_filter_separate_runs(
            sample_ft,
            sample_cols,
            blank_ft,
            blank_cols,
            ratio=args.blank_ratio,
            ppm=args.match_ppm,
            rt_tol=args.match_rt,
        )
    else:
        blank_cols = [c for c in sample_cols if c in blank_names]
        sample_cols = [c for c in sample_cols if c not in blank_names]
        if blank_cols:
            ft, bsum = blanks.blank_filter_same_run(
                sample_ft, sample_cols, blank_cols, args.blank_ratio
            )
            blank_ft = sample_ft[[c for c in sample_ft.columns if c not in sample_cols]]
            ft = ft.drop(columns=blank_cols)
        else:
            ft, bsum, blank_ft = sample_ft, None, None

    if bsum is not None:
        blank_path = os.path.join(tables_dir, "blank_feature_table.tsv")
        blank_ft.to_csv(blank_path, sep="\t", index=False)
        record.add_artifact(
            "blank_feature_table",
            blank_path,
            "Features detected in blank samples",
            "blank_filter",
            n_features=len(blank_ft),
            n_samples=len(blank_cols),
            sample_columns=list(blank_cols),
        )
        record.finish_step("blank_filter", COMPLETED, summary=bsum)
        print(
            f"[pcpfm] Blank filter: {bsum['n_flagged']} of {bsum['n_features']} features flagged "
            f"(blank x {args.blank_ratio} >= sample mean)."
        )
    else:
        record.finish_step(
            "blank_filter", SKIPPED, message="No blank samples identified."
        )

    # feature columns first, then blank columns, then sample intensities in metadata order
    order = {x["sample"]: i for i, x in enumerate(available)}
    sample_cols = sorted(sample_cols, key=lambda c: order.get(c, len(order)))
    lead = [c for c in ft.columns if c not in sample_cols]
    ft = ft[lead + sample_cols]
    ft_path = os.path.join(tables_dir, "feature_table.tsv")
    ft.to_csv(ft_path, sep="\t", index=False)
    record.add_artifact(
        "feature_table",
        ft_path,
        f"Asari {args.asari_table} feature table of non-blank samples, with blank-filter columns",
        "process",
        n_features=len(ft),
        n_samples=len(sample_cols),
        sample_columns=list(sample_cols),
    )
    record.save()
    return {
        "n_inputs": len(inputs),
        "n_files": len(available),
        "n_missing": len(missing),
        "asari_runs": list(runs),
        "n_features": len(ft),
        "n_samples": len(sample_cols),
        "blank_filter": bsum,
    }
