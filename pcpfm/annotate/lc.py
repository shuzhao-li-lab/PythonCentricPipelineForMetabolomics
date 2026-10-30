"""
``pcpfm annotate lc``: annotate an LC-MS feature table.

Levels follow the PCPFM v1 scheme: pre-annotation (empirical compounds) and L4,
L1b (authentic standards, m/z + RT), L2 (MS2 vs a spectral library) and L1a
(authentic standards, m/z + RT + MS2). Each level runs only when its reference
is given. The levels overlap, so their counts are never summed.

Each run writes ``annotations/lc/<run_id>/annotations.tsv`` (long evidence table)
and ``manifest.json`` (provenance). The ``annotations_lc`` artifact points to the
newest successful run; a failed run leaves it unchanged.
"""

import os
import secrets
import time
from argparse import Namespace
from typing import Any

from .. import __version__
from ..project import COMPLETED, FAILED, SKIPPED, ProjectRecord, now
from ..tables import load_table
from .lc_evidence import empty_evidence, file_sha256, write_evidence, write_manifest

WORKFLOW = "lc"
ASARI_WORKFLOW = "LC"
LEVELS = ("grouping", "L4", "L1b", "L2", "L1a")


def source_context(record: ProjectRecord) -> tuple[dict[str, Any], list[str]]:
    """Read the Asari workflow and polarity of the project's sample run.

    Parameters
    ----------
    record : ProjectRecord
        Project record with the processing runs.

    Returns
    -------
    tuple[dict[str, Any], list[str]]
        ``{"workflow", "mode"}`` (either may be None) and warnings.

    Raises
    ------
    ValueError
        If the recorded workflow is not Asari's LC workflow.
    """
    run = record.data.get("runs", {}).get("samples") or {}
    workflow = run.get("workflow")
    # runs recorded before pcpfm kept the Asari mode fall back to the process option
    mode = run.get("mode") or record.step("process").get("params", {}).get("mode")
    warnings = []
    if workflow is None:
        warnings.append(
            "Source Asari workflow is not recorded; assuming LC. "
            "Rerun `pcpfm process` to record it."
        )
    elif workflow != ASARI_WORKFLOW:
        raise ValueError(
            f"The source table comes from the Asari {workflow} workflow; "
            f"`annotate {WORKFLOW}` needs {ASARI_WORKFLOW}."
        )
    return {"workflow": workflow, "mode": mode}, warnings


def run(args: Namespace, record: ProjectRecord) -> None:
    """Run LC annotation and record it under the ``annotate`` step.

    The step summary keeps one entry per workflow, so an LC run does not erase
    the result of another annotation workflow. The step is completed when at
    least one level ran, and skipped otherwise.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for ``annotate lc``.
    record : ProjectRecord
        Project record to update.

    Raises
    ------
    Exception
        Re-raises errors after printing them and marking the step as failed.
    """
    params = {k: v for k, v in vars(args).items() if k not in ("func", "argv")}
    previous = record.step("annotate").get("summary", {}).get("workflows", {})
    record.start_step("annotate", params)
    try:
        summary, warnings = _annotate(args, record, params)
    except Exception as e:
        print(f"[pcpfm] Error: {e}")
        # keep the last successful results of every workflow visible
        record.finish_step(
            "annotate", FAILED, summary={"workflows": previous}, message=str(e)
        )
        raise
    ran = any(lv["status"] == COMPLETED for lv in summary["levels"].values())
    messages = warnings if ran else warnings + ["No annotation level ran."]
    record.finish_step(
        "annotate",
        COMPLETED if ran else SKIPPED,
        summary={"workflows": {**previous, WORKFLOW: summary}},
        message="; ".join(messages) or None,
    )


def _annotate(
    args: Namespace, record: ProjectRecord, params: dict[str, Any]
) -> tuple[dict[str, Any], list[str]]:
    """Validate inputs, run the requested levels and write the run outputs.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed options for ``annotate lc``.
    record : ProjectRecord
        Project record to read and update.
    params : dict[str, Any]
        Resolved options, recorded in the manifest.

    Returns
    -------
    tuple[dict[str, Any], list[str]]
        Workflow summary for the step record, and warnings.
    """
    context, warnings = source_context(record)
    for w in warnings:
        print(f"[pcpfm] Warning: {w}")
    table, sample_cols = load_table(record, args.table)
    source_path = record.artifact_path(args.table)

    run_id = time.strftime("%Y%m%d_%H%M%S") + "_" + secrets.token_hex(4)
    out_dir = record.subdir("annotations", WORKFLOW, run_id)
    levels = {name: {"status": SKIPPED, "reason": "not implemented yet"} for name in LEVELS}
    evidence = empty_evidence()

    evidence_path = os.path.join(out_dir, "annotations.tsv")
    n_rows = write_evidence(evidence, evidence_path)
    source_sha256 = file_sha256(source_path)
    manifest_path = os.path.join(out_dir, "manifest.json")
    write_manifest(
        {
            "workflow": WORKFLOW,
            "run_id": run_id,
            "created": now(),
            "pcpfm_version": __version__,
            "params": params,
            "source": {
                "artifact": args.table,
                "path": record.relpath(source_path),
                "sha256": source_sha256,
                "n_features": len(table),
                "sample_columns": sample_cols,
            },
            "source_context": context,
            "levels": levels,
            "n_evidence_rows": n_rows,
            "outputs": {"evidence": "annotations.tsv"},
            "warnings": warnings,
        },
        manifest_path,
    )
    record.add_artifact(
        "annotations_lc",
        evidence_path,
        "LC annotation evidence (one row per candidate) of the newest successful run",
        "annotate",
        run_id=run_id,
        manifest=record.relpath(manifest_path),
        source_table=args.table,
        source_sha256=source_sha256,
        n_rows=n_rows,
    )
    record.save()
    print(f"[pcpfm] Annotation (lc) run {run_id}: {n_rows} evidence rows in {out_dir}")
    return {"run_id": run_id, "finished": now(), "levels": levels}, warnings
